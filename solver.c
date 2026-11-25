/*
 * solver.c - Cache Lab Bonus：逆向神秘缓存。
 *
 * 求解顺序：
 *   1. b：2 的幂步进，第一次跨出块的偏移就是块大小
 *   2. E 和 RANDOM：步长 2^(b+9) 的同组块装满后探测倒数第二块，
 *      RANDOM 下命中率是 1-1/E（会波动），LRU/FIFO 恒定
 *   3. s：步长 2^(b+d) 装 E+1 个块后回读块 0，首次被逐出的 d 就是 s
 *   4. 写分配：store 一次再 load 同一地址
 *   5. LRU/FIFO：填满一组后触碰块 0 再装新块，回读块 0 区分
 */
#include "solver.h"

#define K 20           /* 比最大的 E=16 大就行 */
#define ROUNDS 120     /* 判断 RANDOM 的轮数 */
#define EST_ROUNDS 4000 /* 估计 RANDOM 的 E 用的轮数 */

void solve(void)
{
    int b, s, E, policy, wa;
    int i, r, d;
    unsigned base_stride, D, S;
    int first = -1, varied = 0, total = 0;

    policy = POLICY_LRU;
    wa = 1;

    /* 1. 块偏移位数 b */
    cache_flush();
    cache_load(0);
    b = 0;
    while (b < 31 && cache_load(1u << b))
        b++;

    /* 2. E 和 RANDOM 判断。
     * 2^(b+9) 是 2^(b+s) 的倍数，i * 2^(b+9) 的所有块都在同一组。装满
     * K 个块后探测倒数第二块：RANDOM 下它只要没被最后一次装入逐出就命中，
     * 命中率是 1-1/E（会波动）；LRU/FIFO 保留最新的 E 块，E>=2 恒命中、
     * E=1 恒缺失。波动则用多轮命中率比对 1-1/E 求 E，恒定则倒序回读
     * 一次直接数出 E。 */
    base_stride = 1u << (b + 9);
    total = 0;
    first = -1;
    varied = 0;
    for (r = 0; r < ROUNDS; r++) {
        int hit;
        cache_flush();
        for (i = 0; i < K; i++)
            cache_load((unsigned)i * base_stride);
        hit = cache_load((unsigned)(K - 2) * base_stride);
        total += hit;
        if (first < 0)
            first = hit;
        else if (hit != first)
            varied = 1;
    }
    if (varied) {
        static const int allowed[9] = {1, 2, 3, 4, 5, 6, 8, 12, 16};
        int best;
        long bestdist, dist, lhs, rhs;
        policy = POLICY_RANDOM;
        for (r = 0; r < EST_ROUNDS - ROUNDS; r++) {
            cache_flush();
            for (i = 0; i < K; i++)
                cache_load((unsigned)i * base_stride);
            total += cache_load((unsigned)(K - 2) * base_stride);
        }
        /* 把实测命中率 total/EST_ROUNDS 和各 E 的理论值 (E-1)/E 比对，
         * 整数化：|total*E - EST_ROUNDS*(E-1)| 最小的就是 E */
        best = allowed[0];
        bestdist = -1;
        for (i = 0; i < 9; i++) {
            lhs = (long)total * allowed[i];
            rhs = (long)EST_ROUNDS * (allowed[i] - 1);
            dist = lhs - rhs;
            if (dist < 0)
                dist = -dist;
            if (bestdist < 0 || dist < bestdist) {
                bestdist = dist;
                best = allowed[i];
            }
        }
        E = best;
    } else {
        /* LRU/FIFO：倒序回读，恰好命中最新装入的 E 块，之后的缺失
         * 只会逐出已经回读过的块，不影响计数 */
        cache_flush();
        for (i = 0; i < K; i++)
            cache_load((unsigned)i * base_stride);
        E = 0;
        for (i = K - 1; i >= 0; i--)
            if (cache_load((unsigned)i * base_stride))
                E++;
    }

    /* 3. 组索引位数 s */
    s = 9; /* 兜底，实际不会走到 */
    for (d = 1; d <= 9; d++) {
        int miss_seen = 0;
        int rounds = (policy == POLICY_RANDOM) ? 12 : 1;
        D = 1u << (b + d);
        for (r = 0; r < rounds; r++) {
            cache_flush();
            for (i = 0; i <= E; i++)
                cache_load((unsigned)i * D);
            if (!cache_load(0))
                miss_seen = 1;
        }
        if (miss_seen) {
            s = d;
            break;
        }
    }

    /* 4. 写分配 */
    cache_flush();
    cache_store(0);
    wa = cache_load(0) ? 1 : 0;

    /* 5. LRU 还是 FIFO（E=1 不计分） */
    if (E > 1 && policy != POLICY_RANDOM) {
        S = 1u << (b + s);
        cache_flush();
        for (i = 0; i < E; i++)
            cache_load((unsigned)i * S);
        cache_load(0);               /* 命中，LRU 下块 0 变最新 */
        cache_load((unsigned)E * S); /* 缺失，逐出一个 */
        policy = cache_load(0) ? POLICY_LRU : POLICY_FIFO;
    }

    report(b, s, E, policy, wa);
}

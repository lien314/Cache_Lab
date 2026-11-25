/*
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 *
 * 分块策略（缓存 s=5, E=1, b=5，块 32B = 8 个 int）：
 *   32x32：8x8 分块，A 的一行先读进 8 个局部变量再写 B，避开对角线冲突
 *   64x64：8x8 块拆成 4x4 象限；对角块借相邻列带做 4x8 暂存
 *   61x67：8x23 分块，边界用带判断的普通循环
 *   其他尺寸：16x16 分块
 */
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/*
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded.
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    int i, bi, bj, ii, jj, k, scr;
    int t0, t1, t2, t3, t4, t5, t6, t7;

    if (M == 32 && N == 32) {
        /* 8x8 分块，A 的一行先读进寄存器再写 B 的列 */
        for (ii = 0; ii < 32; ii += 8) {
            for (jj = 0; jj < 32; jj += 8) {
                for (i = ii; i < ii + 8; i++) {
                    t0 = A[i][jj];     t1 = A[i][jj + 1];
                    t2 = A[i][jj + 2]; t3 = A[i][jj + 3];
                    t4 = A[i][jj + 4]; t5 = A[i][jj + 5];
                    t6 = A[i][jj + 6]; t7 = A[i][jj + 7];
                    B[jj][i] = t0; B[jj + 1][i] = t1;
                    B[jj + 2][i] = t2; B[jj + 3][i] = t3;
                    B[jj + 4][i] = t4; B[jj + 5][i] = t5;
                    B[jj + 6][i] = t6; B[jj + 7][i] = t7;
                }
            }
        }
    } else if (M == 64 && N == 64) {
        for (ii = 0; ii < 64; ii += 8) {
            /* 对角块：先处理，暂存区用相邻列带 */
            scr = (ii == 0) ? 8 : 0;

            /* A 的下半 4x8 搬到暂存区 */
            for (bi = ii; bi < ii + 4; bi++) {
                t0 = A[bi + 4][ii];     t1 = A[bi + 4][ii + 1];
                t2 = A[bi + 4][ii + 2]; t3 = A[bi + 4][ii + 3];
                t4 = A[bi + 4][ii + 4]; t5 = A[bi + 4][ii + 5];
                t6 = A[bi + 4][ii + 6]; t7 = A[bi + 4][ii + 7];
                B[bi][scr] = t0; B[bi][scr + 1] = t1;
                B[bi][scr + 2] = t2; B[bi][scr + 3] = t3;
                B[bi][scr + 4] = t4; B[bi][scr + 5] = t5;
                B[bi][scr + 6] = t6; B[bi][scr + 7] = t7;
            }
            /* 暂存区里两个 4x4 原地转置 */
            for (bi = 0; bi < 4; bi++) {
                for (bj = bi + 1; bj < 4; bj++) {
                    t0 = B[ii + bi][scr + bj];
                    B[ii + bi][scr + bj] = B[ii + bj][scr + bi];
                    B[ii + bj][scr + bi] = t0;
                    t0 = B[ii + bi][scr + bj + 4];
                    B[ii + bi][scr + bj + 4] = B[ii + bj][scr + bi + 4];
                    B[ii + bj][scr + bi + 4] = t0;
                }
            }
            /* A 的上半 4x8 搬到块自身 */
            for (bi = ii; bi < ii + 4; bi++) {
                t0 = A[bi][ii];     t1 = A[bi][ii + 1];
                t2 = A[bi][ii + 2]; t3 = A[bi][ii + 3];
                t4 = A[bi][ii + 4]; t5 = A[bi][ii + 5];
                t6 = A[bi][ii + 6]; t7 = A[bi][ii + 7];
                B[bi][ii] = t0; B[bi][ii + 1] = t1;
                B[bi][ii + 2] = t2; B[bi][ii + 3] = t3;
                B[bi][ii + 4] = t4; B[bi][ii + 5] = t5;
                B[bi][ii + 6] = t6; B[bi][ii + 7] = t7;
            }
            /* 块上半两个 4x4 原地转置 */
            for (bi = ii; bi < ii + 4; bi++) {
                for (bj = bi + 1; bj < ii + 4; bj++) {
                    t0 = B[bi][bj];
                    B[bi][bj] = B[bj][bi];
                    B[bj][bi] = t0;
                    t0 = B[bi][bj + 4];
                    B[bi][bj + 4] = B[bj][bi + 4];
                    B[bj][bi + 4] = t0;
                }
            }
            /* 交换块右上和暂存区左半 */
            for (bi = 0; bi < 4; bi++) {
                t0 = B[ii + bi][ii + 4]; t1 = B[ii + bi][ii + 5];
                t2 = B[ii + bi][ii + 6]; t3 = B[ii + bi][ii + 7];
                B[ii + bi][ii + 4] = B[ii + bi][scr];
                B[ii + bi][ii + 5] = B[ii + bi][scr + 1];
                B[ii + bi][ii + 6] = B[ii + bi][scr + 2];
                B[ii + bi][ii + 7] = B[ii + bi][scr + 3];
                B[ii + bi][scr] = t0; B[ii + bi][scr + 1] = t1;
                B[ii + bi][scr + 2] = t2; B[ii + bi][scr + 3] = t3;
            }
            /* 暂存区整体下移进块下半 */
            for (bi = 0; bi < 4; bi++) {
                B[ii + bi + 4][ii] = B[ii + bi][scr];
                B[ii + bi + 4][ii + 1] = B[ii + bi][scr + 1];
                B[ii + bi + 4][ii + 2] = B[ii + bi][scr + 2];
                B[ii + bi + 4][ii + 3] = B[ii + bi][scr + 3];
                B[ii + bi + 4][ii + 4] = B[ii + bi][scr + 4];
                B[ii + bi + 4][ii + 5] = B[ii + bi][scr + 5];
                B[ii + bi + 4][ii + 6] = B[ii + bi][scr + 6];
                B[ii + bi + 4][ii + 7] = B[ii + bi][scr + 7];
            }

            /* 非对角块 */
            for (jj = 0; jj < 64; jj += 8) {
                if (jj == ii)
                    continue;

                /* 上半 4x8：左半转置到 B 左上，右半暂存到 B 右上 */
                for (bi = jj; bi < jj + 4; bi++) {
                    t0 = A[bi][ii];     t1 = A[bi][ii + 1];
                    t2 = A[bi][ii + 2]; t3 = A[bi][ii + 3];
                    t4 = A[bi][ii + 4]; t5 = A[bi][ii + 5];
                    t6 = A[bi][ii + 6]; t7 = A[bi][ii + 7];
                    B[ii][bi] = t0; B[ii + 1][bi] = t1;
                    B[ii + 2][bi] = t2; B[ii + 3][bi] = t3;
                    B[ii][bi + 4] = t4; B[ii + 1][bi + 4] = t5;
                    B[ii + 2][bi + 4] = t6; B[ii + 3][bi + 4] = t7;
                }
                /* 左下的一列填 B 右上，暂存行下移到 B 左下 */
                for (bj = ii; bj < ii + 4; bj++) {
                    t0 = A[jj + 4][bj]; t1 = A[jj + 5][bj];
                    t2 = A[jj + 6][bj]; t3 = A[jj + 7][bj];
                    t4 = B[bj][jj + 4]; t5 = B[bj][jj + 5];
                    t6 = B[bj][jj + 6]; t7 = B[bj][jj + 7];
                    B[bj][jj + 4] = t0; B[bj][jj + 5] = t1;
                    B[bj][jj + 6] = t2; B[bj][jj + 7] = t3;
                    B[bj + 4][jj] = t4; B[bj + 4][jj + 1] = t5;
                    B[bj + 4][jj + 2] = t6; B[bj + 4][jj + 3] = t7;
                }
                /* 右下象限 */
                for (bi = jj + 4; bi < jj + 8; bi++) {
                    t4 = A[bi][ii + 4]; t5 = A[bi][ii + 5];
                    t6 = A[bi][ii + 6]; t7 = A[bi][ii + 7];
                    B[ii + 4][bi] = t4; B[ii + 5][bi] = t5;
                    B[ii + 6][bi] = t6; B[ii + 7][bi] = t7;
                }
            }
        }
    } else if (M == 61 && N == 67) {
        /* 8x23 分块，边界用带判断的普通循环 */
        for (ii = 0; ii < 67; ii += 8) {
            for (jj = 0; jj < 61; jj += 23) {
                if (ii + 8 <= 67 && jj + 23 <= 61) {
                    for (k = jj; k < jj + 23; k++) {
                        t0 = A[ii][k];     t1 = A[ii + 1][k];
                        t2 = A[ii + 2][k]; t3 = A[ii + 3][k];
                        t4 = A[ii + 4][k]; t5 = A[ii + 5][k];
                        t6 = A[ii + 6][k]; t7 = A[ii + 7][k];
                        B[k][ii] = t0; B[k][ii + 1] = t1;
                        B[k][ii + 2] = t2; B[k][ii + 3] = t3;
                        B[k][ii + 4] = t4; B[k][ii + 5] = t5;
                        B[k][ii + 6] = t6; B[k][ii + 7] = t7;
                    }
                } else {
                    for (i = ii; i < 67 && i < ii + 8; i++) {
                        for (k = jj; k < 61 && k < jj + 23; k++) {
                            B[k][i] = A[i][k];
                        }
                    }
                }
            }
        }
    } else {
        /* 其他尺寸统一用 16x16 分块 */
        for (ii = 0; ii < N; ii += 16) {
            for (jj = 0; jj < M; jj += 16) {
                for (i = ii; i < N && i < ii + 16; i++) {
                    for (k = jj; k < M && k < jj + 16; k++) {
                        B[k][i] = A[i][k];
                    }
                }
            }
        }
    }
}

/*
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }
}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc);

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc);
}

/*
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

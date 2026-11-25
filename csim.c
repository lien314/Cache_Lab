/*
 * csim.c - Cache Lab Part A.
 *
 * Complete access_cache() so the simulator matches the supplied csim-ref.
 *     ./csim [-hv] -s <s> -E <E> -b <b> -t <tracefile>
 *
 * The cache has S = 2^s sets, E lines per set and B = 2^b byte blocks.
 * Replacement is LRU. Assume each access stays within one block.
 * The template parses arguments, allocates and frees the cache, and reads
 * the trace. Replay skips instruction loads (I) and treats a modify (M)
 * as two accesses to the same address. printSummary() prints the totals
 * and writes them to .csim_results.
 */

#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "cachelab.h"

typedef struct {
    int valid;
    uint64_t tag;
    uint64_t last_used;   /* larger means more recently used */
} line_t;

static int s = -1, E = -1, b = -1, verbose = 0;
/* Set i occupies E consecutive lines starting at cache[(size_t)i * E]. */
static line_t *cache = NULL;
static uint64_t tick = 0;
static int hits = 0, misses = 0, evictions = 0;

/* usage - print the help message */
static void usage(const char *prog)
{
    printf("Usage: %s [-hv] -s <num> -E <num> -b <num> -t <file>\n", prog);
    printf("Options:\n");
    printf("  -h         Print this help message.\n");
    printf("  -v         Optional verbose flag.\n");
    printf("  -s <num>   Number of set index bits.\n");
    printf("  -E <num>   Number of lines per set.\n");
    printf("  -b <num>   Number of block offset bits.\n");
    printf("  -t <file>  Trace file.\n");
    printf("\nExamples:\n");
    printf("  linux>  %s -s 4 -E 1 -b 4 -t traces/yi.trace\n", prog);
    printf("  linux>  %s -v -s 8 -E 2 -b 4 -t traces/yi.trace\n", prog);
}

/*
 * 模拟一次对 addr 的数据访问：拆出组索引和标记，按 LRU 更新缓存并
 * 统计 hit / miss / eviction。注意 b 和 s+b 可能等于 64，移位前要判断。
 */
static void access_cache(uint64_t addr)
{
    unsigned i, set;
    uint64_t tag;
    line_t *lines;

    /* set 索引在中间 s 位，标记在高位；b 或 s+b 为 64 时不能直接移位 */
    if (b >= 64) {
        set = 0;
        tag = 0;
    } else {
        set = (unsigned)((addr >> b) & (((uint64_t)1 << s) - 1));
        tag = (s + b >= 64) ? 0 : (addr >> (s + b));
    }

    lines = &cache[(size_t)set * (size_t)E];

    for (i = 0; i < (unsigned)E; i++) {
        if (lines[i].valid && lines[i].tag == tag) {
            hits++;
            lines[i].last_used = tick++;
            if (verbose)
                printf(" hit");
            return;
        }
    }

    misses++;

    /* 先用编号最小的无效行；没有的话换出 last_used 最小的行 */
    unsigned victim = 0;
    int evict = 1;
    for (i = 0; i < (unsigned)E; i++) {
        if (!lines[i].valid) {
            victim = i;
            evict = 0;
            break;
        }
    }
    if (evict) {
        for (i = 1; i < (unsigned)E; i++)
            if (lines[i].last_used < lines[victim].last_used)
                victim = i;
        evictions++;
    }

    lines[victim].valid = 1;
    lines[victim].tag = tag;
    lines[victim].last_used = tick++;

    if (verbose)
        printf(evict ? " miss eviction" : " miss");
}

/* replay - run every data access in the trace through the cache */
static void replay(FILE *fp)
{
    char line[256];
    char op;
    unsigned long long addr;
    int size;

    while (fgets(line, sizeof(line), fp) != NULL) {
        if (sscanf(line, " %c %llx,%d", &op, &addr, &size) != 3)
            continue;
        if (op != 'L' && op != 'S' && op != 'M')
            continue;
        if (verbose)
            printf("%c %llx,%d", op, addr, size);
        access_cache((uint64_t)addr);
        if (op == 'M')
            access_cache((uint64_t)addr);
        if (verbose)
            printf("\n");
    }
}

static int parse_number(const char *text, int option, int minimum)
{
    char *end;
    errno = 0;
    long value = strtol(text, &end, 10);
    if (errno == ERANGE || end == text || *end != '\0' ||
        value < minimum || value > INT_MAX) {
        fprintf(stderr, "Invalid numeric argument for -%c: %s\n", option, text);
        exit(EXIT_FAILURE);
    }
    return (int)value;
}

int main(int argc, char *argv[])
{
    int opt;
    char *tracefile = NULL;
    FILE *fp;

    while ((opt = getopt(argc, argv, "hvs:E:b:t:")) != -1) {
        switch (opt) {
        case 'h':
            usage(argv[0]);
            return 0;
        case 'v':
            verbose = 1;
            break;
        case 's':
            s = parse_number(optarg, opt, 0);
            break;
        case 'E':
            E = parse_number(optarg, opt, 1);
            break;
        case 'b':
            b = parse_number(optarg, opt, 0);
            break;
        case 't':
            tracefile = optarg;
            break;
        default:
            usage(argv[0]);
            return 1;
        }
    }

    if (s < 0 || E <= 0 || b < 0 || tracefile == NULL) {
        printf("%s: Missing required command line argument\n", argv[0]);
        usage(argv[0]);
        return 1;
    }
    if (s > 30 || b > 64 - s) {
        printf("%s: Cache parameters out of range\n", argv[0]);
        return 1;
    }

    cache = calloc((size_t)1 << s, (size_t)E * sizeof(line_t));
    if (cache == NULL) {
        perror("calloc");
        return 1;
    }

    fp = fopen(tracefile, "r");
    if (fp == NULL) {
        perror(tracefile);
        free(cache);
        return 1;
    }
    replay(fp);
    int read_failed = ferror(fp);
    if (read_failed)
        perror(tracefile);
    fclose(fp);
    free(cache);
    if (read_failed)
        return 1;

    printSummary(hits, misses, evictions);
    return 0;
}

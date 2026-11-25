/*
 * mul.c - Matrix multiply C = A * B
 *
 * Each multiply function must have a prototype of the form:
 * void mul(int M, int N, int A[N][M], int B[M][N], int C[N][N]);
 *
 * A multiply function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 *
 * 32x32 的分块策略：三个矩阵的行 i 都落在组 4i..4i+3，A 和 C 的同行
 * 永远冲突，B 第 k 行和 C 第 i 行在 i 与 k 同余 8 时冲突。做法是 8x8
 * 分块，B 的一段行读进 8 个寄存器，C 的一段行读进另外 8 个寄存器，
 * 寄存器里累加完再写回。B 的每个块只取一次，之后都在寄存器里，C 读
 * 改写的时候不会和 B 交替碰同一组。
 */
#include <stdio.h>
#include "cachelab.h"

int is_mul(int M, int N, int A[N][M], int B[M][N], int C[N][N]);

/*
 * mul_submit - This is the solution multiply function that you will be
 *     graded on for Part C of the assignment. Do not change the
 *     description string "Multiply submission", as the driver searches
 *     for that string to identify the multiply function to be graded.
 *     C starts out holding garbage, so every element of C must be written.
 */
char mul_submit_desc[] = "Multiply submission";
void mul_submit(int M, int N, int A[N][M], int B[M][N], int C[N][N])
{
    int i0, j0, k0, i, j, k, tmp;
    int b0, b1, b2, b3, b4, b5, b6, b7;
    int c0, c1, c2, c3, c4, c5, c6, c7;

    if (M == 32 && N == 32) {
        /* 先把 C 清零 */
        for (i = 0; i < 32; i++)
            for (j = 0; j < 32; j++)
                C[i][j] = 0;

        /* 8x8 分块。B 的一段行进 b0..b7，C 的一段行进 c0..c7，寄存器
         * 里累加完再写回；B 的块每个只取一次，避免和 C 的读改写互相逐出 */
        for (i0 = 0; i0 < 32; i0 += 8) {
            for (j0 = 0; j0 < 32; j0 += 8) {
                for (k0 = 0; k0 < 32; k0 += 8) {
                    for (k = k0; k < k0 + 8; k++) {
                        b0 = B[k][j0];     b1 = B[k][j0 + 1];
                        b2 = B[k][j0 + 2]; b3 = B[k][j0 + 3];
                        b4 = B[k][j0 + 4]; b5 = B[k][j0 + 5];
                        b6 = B[k][j0 + 6]; b7 = B[k][j0 + 7];
                        for (i = i0; i < i0 + 8; i++) {
                            c0 = C[i][j0];     c1 = C[i][j0 + 1];
                            c2 = C[i][j0 + 2]; c3 = C[i][j0 + 3];
                            c4 = C[i][j0 + 4]; c5 = C[i][j0 + 5];
                            c6 = C[i][j0 + 6]; c7 = C[i][j0 + 7];
                            c0 += A[i][k] * b0; c1 += A[i][k] * b1;
                            c2 += A[i][k] * b2; c3 += A[i][k] * b3;
                            c4 += A[i][k] * b4; c5 += A[i][k] * b5;
                            c6 += A[i][k] * b6; c7 += A[i][k] * b7;
                            C[i][j0] = c0; C[i][j0 + 1] = c1;
                            C[i][j0 + 2] = c2; C[i][j0 + 3] = c3;
                            C[i][j0 + 4] = c4; C[i][j0 + 5] = c5;
                            C[i][j0 + 6] = c6; C[i][j0 + 7] = c7;
                        }
                    }
                }
            }
        }
    } else {
        /* 其他尺寸的兜底实现 */
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                tmp = 0;
                for (k = 0; k < M; k++)
                    tmp += A[i][k] * B[k][j];
                C[i][j] = tmp;
            }
        }
    }
}

/*
 * mul - A simple baseline multiply function, not optimized for the cache.
 */
char mul_desc[] = "Simple i-j-k multiply";
void mul(int M, int N, int A[N][M], int B[M][N], int C[N][N])
{
    int i, j, k, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            tmp = 0;
            for (k = 0; k < M; k++)
                tmp += A[i][k] * B[k][j];
            C[i][j] = tmp;
        }
    }
}

/*
 * mul_var_kij - 实验用的 k-i-j 分块版本，留作对比
 */
char mul_var_kij_desc[] = "Blocked k-i-j experiment";
void mul_var_kij(int M, int N, int A[N][M], int B[M][N], int C[N][N])
{
    int i0, j0, k0, i, j, k, aik;

    if (M != 32 || N != 32) {
        int tmp;
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++) {
                tmp = 0;
                for (k = 0; k < M; k++)
                    tmp += A[i][k] * B[k][j];
                C[i][j] = tmp;
            }
        return;
    }

    for (i = 0; i < 32; i++)
        for (j = 0; j < 32; j++)
            C[i][j] = 0;

    for (k0 = 0; k0 < 32; k0 += 8) {
        for (j0 = 0; j0 < 32; j0 += 8) {
            for (i0 = 0; i0 < 32; i0 += 8) {
                for (k = k0; k < k0 + 8; k++) {
                    for (i = i0; i < i0 + 8; i++) {
                        aik = A[i][k];
                        for (j = j0; j < j0 + 8; j++) {
                            C[i][j] += aik * B[k][j];
                        }
                    }
                }
            }
        }
    }
}

/*
 * registerFunctions - This function registers your multiply functions
 *     with the driver. At runtime, the driver will evaluate each of the
 *     registered functions and summarize their performance. This is a
 *     handy way to experiment with different multiply strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerMulFunction(mul_submit, mul_submit_desc);

    /* Register any additional multiply functions */
    registerMulFunction(mul, mul_desc);
    registerMulFunction(mul_var_kij, mul_var_kij_desc);
}

/*
 * is_mul - This helper function checks if C is the product of A and B.
 *     You can check the correctness of your multiply by calling it
 *     before returning from the multiply function.
 */
int is_mul(int M, int N, int A[N][M], int B[M][N], int C[N][N])
{
    int i, j, k, num;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            num = 0;
            for (k = 0; k < M; k++)
                num += A[i][k] * B[k][j];
            if (num != C[i][j])
                return 0;
        }
    }
    return 1;
}

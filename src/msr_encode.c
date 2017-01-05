//
// Created by syd on 16-11-1.
//
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "msr.h"
#include "gf.h"

#include <immintrin.h>


#define MAX_NODE 16
#define MAX_STRIPE 64


#define AVX2

#include <stdint.h>
#include <assert.h>
#include "gf.h"


#ifdef AVX2

typedef __m256i encode_t;

#define REGION_SIZE 512
#define REGION_BLOCKS (REGION_SIZE/sizeof(encode_t))


__m256i mask_lo, mask_hi;

__m256i low_table[256];
__m256i high_table[256];

void init_avx2() {
    mask_lo = _mm256_set1_epi8(0x0f);
    mask_hi = _mm256_set1_epi8(0xf0);
    //print_encode(mask_lo);
    //print_encode(mask_hi);

    int low_array[32];
    int high_array[32];
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 16; j++)
            low_array[j] = gf_mul(i, j);
        for (int j = 0; j < 16; j++)
            low_array[j + 16] = gf_mul(i, j);
        for (int j = 0; j < 16; j++)
            high_array[j] = gf_mul(i, j << 4);
        for (int j = 0; j < 16; j++)
            high_array[j + 16] = gf_mul(i, j << 4);
        low_table[i] = _mm256_setr_epi8(low_array[0], low_array[1], low_array[2], low_array[3], low_array[4],
                                        low_array[5], low_array[6], low_array[7], low_array[8], low_array[9],
                                        low_array[10], low_array[11], low_array[12], low_array[13], low_array[14],
                                        low_array[15], low_array[16], low_array[17], low_array[18], low_array[19],
                                        low_array[20], low_array[21], low_array[22], low_array[23], low_array[24],
                                        low_array[25], low_array[26], low_array[27], low_array[28], low_array[29],
                                        low_array[30], low_array[31]);
        high_table[i] = _mm256_setr_epi8(high_array[0], high_array[1], high_array[2], high_array[3], high_array[4],
                                         high_array[5], high_array[6], high_array[7], high_array[8], high_array[9],
                                         high_array[10], high_array[11], high_array[12], high_array[13], high_array[14],
                                         high_array[15], high_array[16], high_array[17], high_array[18], high_array[19],
                                         high_array[20], high_array[21], high_array[22], high_array[23], high_array[24],
                                         high_array[25], high_array[26], high_array[27], high_array[28], high_array[29],
                                         high_array[30], high_array[31]);

        //print_encode(low_table[i]);
        //print_encode(high_table[i]);
    }


}

inline __m256i xor_region(__m256i input1, __m256i input2) {

    return _mm256_xor_si256(input1, input2);


}

inline __m256i multiply_region(__m256i input, uint8_t x) {

    __m256i low = _mm256_and_si256(input, mask_lo);
    __m256i high = _mm256_and_si256(input, mask_hi);

    __m256i low_t = low_table[x], high_t = high_table[x];
    high = _mm256_srli_epi16(high, 4);

    __m256i right = _mm256_shuffle_epi8(low_t, low);
    __m256i left = _mm256_shuffle_epi8(high_t, high);


    return _mm256_xor_si256(left, right);

}

void print_encode(__m256i x) {
    //for(int i=0;i<32;i++)
    printf("%0x ", _mm256_extract_epi8(x, 0));
    printf("%0x ", _mm256_extract_epi8(x, 1));
    printf("%0x ", _mm256_extract_epi8(x, 2));
    printf("%0x ", _mm256_extract_epi8(x, 3));
    printf("%0x ", _mm256_extract_epi8(x, 4));
    printf("%0x ", _mm256_extract_epi8(x, 5));
    printf("%0x ", _mm256_extract_epi8(x, 6));
    printf("%0x ", _mm256_extract_epi8(x, 7));
    printf("%0x ", _mm256_extract_epi8(x, 8));
    printf("%0x ", _mm256_extract_epi8(x, 9));
    printf("%0x ", _mm256_extract_epi8(x, 10));

    printf("\n");
}


#endif


uint32_t GfPow[GF_SIZE << 1];
uint32_t GfLog[GF_SIZE << 1];

void gf_init() {
    uint32_t i;

    GfPow[0] = 1;
    GfLog[0] = GF_SIZE;
    for (i = 1; i < GF_SIZE; i++) {
        GfPow[i] = GfPow[i - 1] << 1;
        if (GfPow[i] >= GF_SIZE) {
            GfPow[i] ^= GF_MOD;
        }
        GfPow[i + GF_SIZE - 1] = GfPow[i];
        GfLog[GfPow[i] + GF_SIZE - 1] = GfLog[GfPow[i]] = i;
    }
    init_avx2();
}

inline uint8_t gf_mul(uint32_t a, uint32_t b) {
    if (a && b) {
        return GfPow[GfLog[a] + GfLog[b]];
    }

    return 0;
}

inline uint8_t gf_div(uint32_t a, uint32_t b) {
    if (b) {
        if (a) {
            return GfPow[GF_SIZE - 1 + GfLog[a] - GfLog[b]];
        }
        return 0;
    }
    assert(0);
}

uint32_t _pow(uint32_t a, int b) {
    uint32_t ret = 1;
    for (int i = 1; i <= b; i++)
        ret *= a;
    return ret;
}


inline int get_bit(int z, int y, int q, int t) {
    return z / _pow(q, t - y - 1) % q;
}


int permute(int z, int remove_bit_id, int added_bit, int q, int t) {
    int result = 0;
    int power = _pow(q, t - 1);
    for (int i = 0; i < t; i++) {
        if (i != remove_bit_id)
            result = result * q + (z / power) % q;
        else
            result = result * q + added_bit;
        power /= q;
    }
    return result;
}

uint8_t gf_pow(int x, int y) {
    if (y == 0)
        return 1;
    if (y == 1)
        return x;
    if (y % 2 == 0) {
        uint8_t mid = gf_pow(x, y / 2);
        return gf_mul(mid, mid);
    } else {
        uint8_t mid = gf_pow(x, y / 2);
        return gf_mul(x, gf_mul(mid, mid));
    }
}

const uint8_t u = 3;

static uint8_t node_companion[MAX_NODE][MAX_STRIPE];
static uint8_t z_companion[MAX_NODE][MAX_STRIPE];
static uint8_t theta[MAX_NODE][MAX_NODE];
static uint8_t u_theta[MAX_NODE][MAX_NODE];

static uint8_t inv_matrix[MAX_NODE][MAX_NODE];
static uint8_t final_matrix[MAX_NODE][MAX_NODE];
static uint8_t u_final_matrix[MAX_NODE][MAX_NODE];
static uint8_t matrix[MAX_NODE][MAX_NODE];

static void invert_matrix_rs(int *errors, int error_cnt, int n, int k) {

    memset(inv_matrix, 0, sizeof(inv_matrix));

    for (int i = 0; i < error_cnt; i++)
        inv_matrix[i][i] = 1;

    bool is_error[n];

    memset(is_error, 0, sizeof(bool) * n);

    for (int i = 0; i < error_cnt; i++)
        is_error[errors[i]] = 1;
    int j;

    memset(matrix, 0, sizeof(matrix));
    for (int i = 0; i < error_cnt; i++)
        for (j = 0; j < error_cnt; j++) {
            matrix[i][j] = theta[i][errors[j]];
        }


    //Gaussian Elimination
    for (int i = 0; i < error_cnt; i++) {

        if (!matrix[i][i]) {
            for (j = i + 1; j < error_cnt; j++)
                if (matrix[j][i])
                    break;
            assert(j != error_cnt);

            for (int t = 0; t < error_cnt; t++) {
                uint8_t tmp = matrix[i][t];
                matrix[i][t] = matrix[j][t];
                matrix[j][t] = tmp;

                tmp = inv_matrix[i][t];
                inv_matrix[i][t] = inv_matrix[j][t];
                inv_matrix[j][t] = tmp;

            }

        }

        uint8_t f = matrix[i][i];

        for (j = 0; j < error_cnt; j++) {
            matrix[i][j] = gf_div(matrix[i][j], f);
            inv_matrix[i][j] = gf_div(inv_matrix[i][j], f);
        }

        for (j = 0; j < error_cnt; j++)
            if (i != j) {
                f = matrix[j][i];
                for (int t = 0; t < error_cnt; t++) {
                    matrix[j][t] ^= gf_mul(matrix[i][t], f);
                    inv_matrix[j][t] ^= gf_mul(inv_matrix[i][t], f);
                }
            }
    }

    memset(final_matrix, 0, sizeof(final_matrix));


    for (int i = 0; i < error_cnt; i++)
        for (j = 0; j < n; j++) {

            for (int t = 0; t < error_cnt; t++) {
                final_matrix[i][j] ^= gf_mul(inv_matrix[i][t], theta[t][j]);
            }

        }

    for (int i = 0; i < error_cnt; i++) {
        for (j = 0; j < n; j++) {
            printf("%u ", final_matrix[i][j]);
        }
        printf("\n");
    }

    for (int i = 0; i < error_cnt; i++)
        for (j = 0; j < n; j++) {
                u_final_matrix[i][j] = gf_mul(final_matrix[i][j], u);
        }


}

static void invert_matrix_regenerate(int x0, int y0, int q, int t) {
    memset(inv_matrix, 0, sizeof(inv_matrix));

    for (int i = 0; i < q; i++)
        inv_matrix[i][i] = 1;

    uint8_t matrix[q][q];
    for (int i = 0; i < q; i++)
        for (int j = 0; j < q; j++) {
            int id = y0 * q + j;
            if (j == x0)
                matrix[i][j] = theta[i][id];
            else
                matrix[i][j] = gf_mul(u, theta[i][id]);
        }

    for (int i = 0; i < q; i++) {
        uint8_t f = matrix[i][i];

        for (int j = 0; j < q; j++) {
            matrix[i][j] = gf_div(matrix[i][j], f);
            inv_matrix[i][j] = gf_div(inv_matrix[i][j], f);
        }

        for (int j = 0; j < q; j++)
            if (i != j) {
                f = matrix[j][i];
                for (int k = 0; k < q; k++) {
                    matrix[j][k] ^= gf_mul(matrix[i][k], f);
                    inv_matrix[j][k] ^= gf_mul(inv_matrix[i][k], f);
                }
            }
    }
}


static void init_companion(int n, int k) {
    assert(n <= MAX_NODE);
    int q = n - k;
    int t = n / q;
    int stripe = _pow(q, t);
    assert(stripe <= MAX_STRIPE);

    for (int i = 0; i < n; i++)
        for (int z = 0; z < stripe; z++) {
            int x = i % q;
            int y = i / q;
            node_companion[i][z] = get_bit(z, y, q, t) + y * q;
            z_companion[i][z] = permute(z, y, x, q, t);
        }

}

static void init_theta(int n, int k) {
    memset(theta, 0, sizeof(theta));

    for (int i = 0; i < n - k; i++)
        for (int j = 0; j < k; j++) {
            theta[i][j] = gf_div(1, j ^ (i + k));
            u_theta[i][j] = gf_mul(u, theta[i][j]);
        }

    for (int i = 0; i < n - k; i++) {
        theta[i][i + k] = 1;
        u_theta[i][i + k] = u;
    }

}

void init(int n, int k) {
    init_companion(n, k);
    init_theta(n, k);
}

int test_companion(int i, int z, int q, int t) {
    int x = i % q;
    int y = i / q;

    int companion = get_bit(z, y, q, t) + y * q;
    int new_z = permute(z, y, x, q, t);


    int com_x = companion % q;
    int com_y = companion / q;

    int com_com = get_bit(new_z, com_y, q, t) + com_y * q;
    assert(com_com == i);
    int com_com_z = permute(new_z, com_y, com_x, q, t);
    assert(com_com_z == z);

}


static void solve_equation(int *errors, int error_cnt, uint8_t *kappa) {


    encode_t res[error_cnt];
    memset(res, 0, sizeof(res));
    uint8_t *res_ptr = (uint8_t *) res;

    for (int j = 0; j < error_cnt; j++) {
        //res[j] = 0;
        for (int k = 0; k < error_cnt; k++)
            for (int w = 0; w < sizeof(encode_t); w++)
                res_ptr[j * sizeof(encode_t) + w] ^= gf_mul(inv_matrix[j][k], kappa[k * sizeof(encode_t) + w]);
    }

    for (int j = 0; j < error_cnt; j++)
        ((encode_t *) kappa)[j] = res[j];

}

//the number of each node range from 0 to n-1.
static int compute_sigma(int *errors, int error_cnt, int q, int t, int z) {
    int sigma = 0;
    for (int i = 0; i < error_cnt; i++) {
        if (errors[i] % q == get_bit(z, errors[i] / q, q, t)) {
            sigma++;
        }
    }
    return sigma;
}

encode_t *data_ptr[MAX_NODE];
encode_t kappa[MAX_NODE][REGION_BLOCKS + 3];

int msr_regenerate(int len, int n, int k, uint8_t **data_collected, uint8_t *memory) {
    int i;
    int q = n - k;
    int t = n / q;
    int stripe_size = _pow(q, t);
    int beta = stripe_size / q;


    int block_size = len / beta / REGION_SIZE;
    int index;

    int error_id = -1;
    for (i = 0; i < n; i++)
        if (!data_collected[i]) {
            error_id = i;
        }

    assert(error_id != -1);

    int y_0 = error_id / q;
    int x_0 = error_id % q;

    encode_t *mem_ptr = (encode_t *) memory;

    for (i = 0; i < n; i++)
        data_ptr[i] = (encode_t *) data_collected[i];

    invert_matrix_regenerate(x_0, y_0, q, t);

    int z_num[MAX_STRIPE];
    int z_pos[MAX_STRIPE];
    memset(z_pos, -1, sizeof(z_pos));
    for (int z_id = 0; z_id < beta; z_id++) {

        z_num[z_id] = (z_id / _pow(q, t - y_0 - 1) * q + x_0) * _pow(q, t - y_0 - 1) + z_id % _pow(q, t - y_0 - 1);
        z_pos[z_num[z_id]] = z_id;
    }

    for (index = 0; index < block_size; index++) {
        for (int z_id = 0; z_id < beta; z_id++) {

            for (int j = 0; j < q; j++)
                for (int w = 0; w < REGION_BLOCKS; w++)
                    kappa[j][w] = _mm256_set1_epi8(0);


            int z = z_num[z_id];

            for (i = 0; i < n; i++)
                if (data_collected[i]) {
                    for (int e = 0; e < q; e++)
                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            kappa[e][w] = xor_region(kappa[e][w], multiply_region(
                                    data_ptr[i][(z_id * block_size + index) * REGION_BLOCKS + w], theta[e][i]));
                        }
                }

            for (i = 0; i < n; i++) {
                int companion = node_companion[i][z];
                if (i != companion && data_collected[companion]) {
                    int z_comp = z_companion[i][z];
                    z_comp = z_pos[z_comp];
                    assert(z_comp != -1);

                    for (int e = 0; e < q; e++)
                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            kappa[e][w] = xor_region(kappa[e][w], multiply_region(
                                    data_ptr[companion][(block_size * z_comp + index) * REGION_BLOCKS + w],
                                    u_theta[e][i]));
                        }
                }
            }


            for (int j = 0; j < q; j++) {

                int z = (z_id / _pow(q, t - y_0 - 1) * q + j) * _pow(q, t - y_0 - 1) + z_id % _pow(q, t - y_0 - 1);
                int z_index = (z * block_size + index) * REGION_BLOCKS;

                for (int w = 0; w < REGION_BLOCKS; w++)
                    mem_ptr[z_index + w] = _mm256_set1_epi8(0);

                for (int k = 0; k < q; k++)
                    for (int w = 0; w < REGION_BLOCKS; w++)
                        mem_ptr[z_index + w] = xor_region(mem_ptr[z_index + w],
                                                          multiply_region(kappa[k][w],
                                                                          inv_matrix[j][k]));
            }

        }
    }
}


static void
sequential_decode(uint8_t **data, int index, int block_size, int *errors, int error_cnt,
                  int *sigmas,
                  int q, int t) {

    //clock_t start = clock();
    //printf("Sequential Decode:\n");
    int n = q * t;
    uint32_t z_total = _pow(q, t);
    int K = n - q;

    int sigma_max = 0;
    for (int z = 0; z < z_total; z++) {

        if (sigmas[z] > sigma_max)
            sigma_max = sigmas[z];
    }

    //printf("Time for compute sigma:%lf\n",(clock() - start)/(double)(CLOCKS_PER_SEC));

    bool is_error[n];

    memset(is_error, 0, sizeof(bool) * n);

    for (int i = 0; i < error_cnt; i++)
        is_error[errors[i]] = 1;

    uint8_t a = gf_div(1, gf_mul(1 ^ u, 1 ^ u));
    uint8_t b = gf_div(u, gf_mul(1 ^ u, 1 ^ u));

    assert(a ^ gf_mul(b, u) == 1);
    assert(b ^ gf_mul(a, u) == u);

    int s = 0;


    for (int i = 0; i < n; i++)
        data_ptr[i] = (encode_t *) data[i];

    while (s <= sigma_max) {
        for (int z = 0; z < z_total; z++)
            if (sigmas[z] == s) {
                //printf("%d\n",z);
                int z_index = (block_size * z + index) * REGION_BLOCKS;


                for (int j = 0; j < error_cnt; j++)
                    for (int w = 0; w < REGION_BLOCKS; w++)
                        kappa[j][w] = _mm256_set1_epi8(0);


                for (int i = 0; i < n; i++) {
                    if (!is_error[i]) {

                        for (int j = 0; j < error_cnt; j++)
                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                kappa[j][w] = xor_region(kappa[j][w],
                                                         multiply_region(data_ptr[i][z_index + w], theta[j][i]));
                            }
                    }
                    int companion = node_companion[i][z];
                    if (i != companion && (!is_error[companion] || !is_error[i])) {
                        int z_comp = z_companion[i][z];

                        for (int j = 0; j < error_cnt; j++)
                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                kappa[j][w] = xor_region(kappa[j][w], multiply_region(
                                        data_ptr[companion][(block_size * z_comp + index) * REGION_BLOCKS + w],
                                        u_theta[j][i]));
                            }
                    }

                }


                for (int j = 0; j < error_cnt; j++) {

                    for (int w = 0; w < REGION_BLOCKS; w++)
                        data_ptr[errors[j]][z_index + w] = _mm256_set1_epi8(0);

                    for (int k = 0; k < error_cnt; k++)
                        for (int w = 0; w < REGION_BLOCKS; w++)
                            data_ptr[errors[j]][z_index + w] = xor_region(data_ptr[errors[j]][z_index + w],
                                                                          multiply_region(kappa[k][w],
                                                                                          inv_matrix[j][k]));
                }

            }


        for (int z = 0; z < z_total; z++)
            if (sigmas[z] == s) {

                for (int j = 0; j < error_cnt; j++) {

                    int error = errors[j];

                    int companion = node_companion[error][z];
                    int new_z = z_companion[error][z];

                    if (companion < error && is_error[companion]) {

                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            int z_pos = (block_size * z + index) * REGION_BLOCKS + w;
                            int new_z_pos = (block_size * new_z + index) * REGION_BLOCKS + w;

                            encode_t a_cur = xor_region(multiply_region(data_ptr[error][z_pos], a),
                                                        multiply_region(data_ptr[companion][new_z_pos], b));

                            encode_t a_companion = xor_region(multiply_region(data_ptr[error][z_pos], b),
                                                              multiply_region(data_ptr[companion][new_z_pos], a));


                            data_ptr[error][z_pos] = a_cur;
                            data_ptr[companion][new_z_pos] = a_companion;
                        }

                    }

                }


            }


        s++;

    }
}

encode_t tmp[MAX_NODE][MAX_STRIPE * REGION_BLOCKS + 3];

int err_id[MAX_NODE];


static void
sequential_decode_fast(uint8_t **data, int index, int block_size, int *errors, int error_cnt,
                       int *sigmas,
                       int q, int t) {
    int n = q * t;
    int k = n - q;

    int stripe_size = _pow(q, t);

    int sigma_max = 0;
    for (int z = 0; z < stripe_size; z++) {
        if (sigmas[z] > sigma_max)
            sigma_max = sigmas[z];
    }

    bool is_error[n];

    memset(is_error, 0, sizeof(bool) * n);

    for (int i = 0; i < error_cnt; i++)
        is_error[errors[i]] = 1;

    uint8_t a = gf_div(1, gf_mul(1 ^ u, 1 ^ u));
    uint8_t b = gf_div(u, gf_mul(1 ^ u, 1 ^ u));

    assert(a ^ gf_mul(b, u) == 1);
    assert(b ^ gf_mul(a, u) == u);

    int s = 0;
    //memcpy(data_buffer,data_chunk,sizeof(data_buffer));

    encode_t *data_ptr[n];
    for (int i = 0; i < n; i++)
        data_ptr[i] = (encode_t *) data[i];

    //memset(res,0,sizeof(res));
    for (int i = 0; i < error_cnt; i++)
        for (int z = 0; z < stripe_size * REGION_BLOCKS; z++)
            tmp[i][z] = _mm256_setzero_si256();

    for (int i = 0; i < error_cnt; i++)
        err_id[errors[i]] = i;

    int ok[k + 1];
    ok[k] = 0;
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (!is_error[i])
            ok[j++] = i;
    }


    //The construction of B.
    while (s <= sigma_max) {
        for (int z = 0; z < stripe_size; z++)
            if (sigmas[z] == s) {

                for (j = 0; j < k; j++) {
                    int node_id = ok[j];

                    int z_index = (block_size * z + index) * REGION_BLOCKS;

                    int companion = node_companion[node_id][z];

                    if (companion != node_id) {
                        int new_z = z_companion[node_id][z];
                        int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            encode_t a_cur = xor_region(data_ptr[node_id][z_index + w],
                                                        multiply_region(data_ptr[companion][new_z_index + w], u));


                            for (int e = 0; e < 4; e++)
                                tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                           multiply_region(a_cur,
                                                                                           final_matrix[e][node_id]));


                            _mm_prefetch(&data_ptr[ok[j + 1]][(block_size * z + index) * REGION_BLOCKS + w],
                                         _MM_HINT_NTA);

                            _mm_prefetch(&data_ptr[node_companion[ok[j + 1]][z]][
                                                 (block_size * z_companion[ok[j + 1]][z] + index) * REGION_BLOCKS + w],
                                         _MM_HINT_NTA);


                        }


                    } else if (companion == node_id) {
                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            for (int e = 0; e < 4; e++) {

                                tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                           multiply_region(data_ptr[node_id][
                                                                                                   (block_size * z +
                                                                                                    index) *
                                                                                                   REGION_BLOCKS + w],
                                                                                           final_matrix[e][node_id]));

                            }
                            _mm_prefetch(&data_ptr[ok[j + 1]][(block_size * z + index) * REGION_BLOCKS + w],
                                         _MM_HINT_NTA);
                            _mm_prefetch(&data_ptr[node_companion[ok[j + 1]][z]][
                                                 (block_size * z_companion[ok[j + 1]][z] + index) * REGION_BLOCKS + w],
                                         _MM_HINT_NTA);

                        }

                    }
                }

                for (j = 0; j < error_cnt; j++) {
                    int error = errors[j];

                    int companion = node_companion[error][z];
                    if (error != companion) {
                        int new_z = z_companion[error][z];
                        int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            for (int e = 0; e < 4; e++) {
                                tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                           multiply_region(data_ptr[companion][
                                                                                                   new_z_index + w],
                                                                                                  u_final_matrix[e][error]));
                            }

                            if(j != error_cnt -1)
                            _mm_prefetch(&data_ptr[node_companion[errors[j + 1]][z]][
                                                 (block_size * z_companion[errors[j + 1]][z] + index) * REGION_BLOCKS + w],
                                         _MM_HINT_NTA);


                        }
                    }

                }


            }

        for (int z = 0; z < stripe_size; z++)
            if (sigmas[z] == s)

                for (j = 0; j < error_cnt; j++) {

                    int error = errors[j];

                    int companion = node_companion[error][z];
                    int new_z = z_companion[error][z];

                    int comp_id = err_id[companion];

                    if (companion < error && is_error[companion]) {
                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            int z_index = z * REGION_BLOCKS + w;
                            int new_z_index = new_z * REGION_BLOCKS + w;

                            encode_t a_cur
                                    = xor_region(multiply_region(tmp[j][z_index], a),
                                                 multiply_region(
                                                         tmp[comp_id][new_z_index], b));

                            encode_t a_companion = xor_region(
                                    multiply_region(tmp[j][z_index], b),
                                    multiply_region(tmp[comp_id][new_z_index], a));


                            _mm256_stream_si256(&data_ptr[error][(block_size * z + index) * REGION_BLOCKS + w], a_cur);
                            _mm256_stream_si256(&data_ptr[companion][(block_size * new_z + index) * REGION_BLOCKS + w],
                                                a_companion);
                        }
                    } else if (companion == error || !is_error[companion]) {
                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            _mm256_stream_si256(&data_ptr[error][(block_size * z + index) * REGION_BLOCKS + w],
                                                tmp[j][z * REGION_BLOCKS + w]);
                        }
                    }
                }

        s++;
    }

}


static void
systematic_encode_with_prefetch(uint8_t **data, int index, int block_size, int *errors, int error_cnt,
                                int *sigmas,
                                int q, int t) {

    int n = q * t;
    int k = n - q;

    int stripe_size = _pow(q, t);

    int sigma_max = 0;
    for (int z = 0; z < stripe_size; z++) {
        if (sigmas[z] > sigma_max)
            sigma_max = sigmas[z];
    }

    bool is_error[n];

    memset(is_error, 0, sizeof(bool) * n);

    for (int i = 0; i < error_cnt; i++)
        is_error[errors[i]] = 1;

    uint8_t a = gf_div(1, gf_mul(1 ^ u, 1 ^ u));
    uint8_t b = gf_div(u, gf_mul(1 ^ u, 1 ^ u));

    assert(a ^ gf_mul(b, u) == 1);
    assert(b ^ gf_mul(a, u) == u);

    int s = 0;
    //memcpy(data_buffer,data_chunk,sizeof(data_buffer));

    encode_t *data_ptr[n];
    for (int i = 0; i < n; i++)
        data_ptr[i] = (encode_t *) data[i];

    //memset(res,0,sizeof(res));
    for (int i = 0; i < error_cnt; i++)
        for (int z = 0; z < stripe_size * REGION_BLOCKS; z++)
            tmp[i][z] = _mm256_set1_epi8(0);

    for (int i = 0; i < error_cnt; i++)
        err_id[errors[i]] = i;

    //The construction of B.
    while (s <= sigma_max) {
        for (int z = 0; z < stripe_size; z++)
            if (sigmas[z] == s) {

                for (int j = 0; j < k; j++) {

                    int z_index = (block_size * z + index) * REGION_BLOCKS;

                    int companion = node_companion[j][z];
                    if (companion != j && !is_error[companion]) {

                        int new_z = z_companion[j][z];


                        int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;


                        for (int w = 0; w < REGION_BLOCKS; w++) {

                            encode_t a_cur = xor_region(data_ptr[j][z_index + w],
                                                        multiply_region(data_ptr[companion][new_z_index + w], u));


                            for (int e = 0; e < 4; e++)
                                tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                           multiply_region(a_cur,
                                                                                           theta[errors[e] - k][j]));


                            _mm_prefetch(&data_ptr[j + 1][(block_size * z + index) * REGION_BLOCKS + w], _MM_HINT_NTA);
                            _mm_prefetch(&data_ptr[node_companion[j + 1][z]][
                                                 (block_size * z_companion[j + 1][z] + index) * REGION_BLOCKS + w],
                                         _MM_HINT_NTA);

                        }

                    } else if (companion == j) {
                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            for (int e = 0; e < 4; e++) {
                                tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                           multiply_region(data_ptr[j][
                                                                                                   (block_size * z +
                                                                                                    index) *
                                                                                                   REGION_BLOCKS + w],
                                                                                           theta[errors[e] - k][j]));
                                //tmp[e][z*REGION_BLOCKS + w] = mid;
                                //tmp[errors[e]][(block_size * z + index) * REGION_BLOCKS + w]=mid;

                                //_mm256_stream_si256(&data_ptr[errors[e]][(block_size * z + index) * REGION_BLOCKS + w],mid);
                            }
                            _mm_prefetch(&data_ptr[j + 1][(block_size * z + index) * REGION_BLOCKS + w], _MM_HINT_NTA);
                            _mm_prefetch(&data_ptr[node_companion[j + 1][z]][
                                                 (block_size * z_companion[j + 1][z] + index) * REGION_BLOCKS + w],
                                         _MM_HINT_NTA);

                        }
                        //print_encode(tmp[0][z*REGION_BLOCKS + w]);

                    }

                }
            }


        for (int z = 0; z < stripe_size; z++)
            if (sigmas[z] == s)

                for (int j = 0; j < error_cnt; j++) {

                    int error = errors[j];

                    int companion = node_companion[error][z];
                    int new_z = z_companion[error][z];

                    int comp_id = err_id[companion];
                    //printf("%d %d\n",companion,error);

                    if (companion < error && is_error[companion]) {
                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            int z_index = z * REGION_BLOCKS + w;
                            int new_z_index = new_z * REGION_BLOCKS + w;

                            encode_t a_cur
                                    = xor_region(multiply_region(tmp[j][z_index], a),
                                                 multiply_region(
                                                         tmp[comp_id][new_z_index], b));

                            encode_t a_companion = xor_region(
                                    multiply_region(tmp[j][z_index], b),
                                    multiply_region(tmp[comp_id][new_z_index], a));


                            _mm256_stream_si256(&data_ptr[error][(block_size * z + index) * REGION_BLOCKS + w], a_cur);
                            _mm256_stream_si256(&data_ptr[companion][(block_size * new_z + index) * REGION_BLOCKS + w],
                                                a_companion);
                        }
                    } else if (companion == error) {
                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            _mm256_stream_si256(&data_ptr[error][(block_size * z + index) * REGION_BLOCKS + w],
                                                tmp[j][z * REGION_BLOCKS + w]);
                        }
                    }
                }
        s++;
    }

}


void msr_encode(int len, int n, int k, uint8_t **data, uint8_t **memory_allocated) {
    int r = n - k;

    assert(n % r == 0);

    int q = r;
    int t = n / q;


    uint32_t stripe_size = _pow(q, t);

    assert(len % (stripe_size * REGION_SIZE) == 0);

    int error_cnt = 0;

    for (int i = 0; i < n; i++)
        if (!data[i])
            error_cnt++;


    int errors[error_cnt];

    int index = 0;
    for (int i = 0; i < n; i++)
        if (!data[i])
            errors[index++] = i;

    bool is_error[n];

    memset(is_error, 0, sizeof(bool) * n);

    bool is_systematic = true;

    for (int i = 0; i < error_cnt; i++) {
        is_error[errors[i]] = 1;
        if (errors[i] < k)
            is_systematic = false;
    }

    index = 0;

    for (int i = 0; i < n; i++)
        if (is_error[i]) {
            data[i] = memory_allocated[index++];
        }

    invert_matrix_rs(errors, error_cnt, n, k);

    //is_systematic = false;
    int sigmas[stripe_size];

    for (int z = 0; z < stripe_size; z++)
        sigmas[z] = compute_sigma(errors, error_cnt, q, t, z);

    int block_size = len / stripe_size / REGION_SIZE;


    for (index = 0; index < block_size; index++) {
        if (is_systematic)
            sequential_decode_fast(data, index, block_size, errors, error_cnt, sigmas, q, t);
        else
            sequential_decode_fast(data, index, block_size, errors, error_cnt, sigmas, q, t);
    }
}
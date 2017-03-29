//
// Created by syd on 16-11-1.
//
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <stdbool.h>
#include "msr.h"
#include "gf.h"

#include <immintrin.h>
#include <time.h>


#define MAX_NODE 16
#define MAX_STRIPE 64


#define AVX2


#ifdef AVX2

typedef __m256i encode_t;


#define REGION_SIZE 512
#define REGION_BLOCKS (REGION_SIZE/sizeof(encode_t))

uint32_t _pow(uint32_t a, int b) {
    uint32_t ret = 1;
    for (int i = 1; i <= b; i++)
        ret *= a;
    return ret;
}


static __m256i mask_lo, mask_hi;

static __m256i low_table[256];
static __m256i high_table[256];

void init_avx2() {
    mask_lo = _mm256_set1_epi8(0x0f);
    mask_hi = _mm256_set1_epi8(0xf0);

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

inline static __m256i multiply_region(__m256i input, uint8_t x) {

    __m256i low = _mm256_and_si256(input, mask_lo);
    __m256i high = _mm256_and_si256(input, mask_hi);

    __m256i low_t = low_table[x], high_t = high_table[x];
    //high = _mm256_srli_epi64(high,4);
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

    printf("\n");
}


#endif


static uint32_t GfPow[GF_SIZE << 1];
static uint32_t GfLog[GF_SIZE << 1];

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

//Read only after the init, so can be placed globally thread-safely.
static uint8_t node_companion[MAX_NODE][MAX_STRIPE];
static uint8_t z_companion[MAX_NODE][MAX_STRIPE];
static uint8_t theta[MAX_NODE][MAX_NODE];
static uint8_t u_theta[MAX_NODE][MAX_NODE];


//For DFS.
int fill[MAX_NODE][MAX_STRIPE];
int pointer[MAX_NODE];
char used[MAX_NODE][MAX_STRIPE];
int reorder[MAX_NODE][MAX_STRIPE];

bool dfs(int i,int z){
    int next_i = (i+1) % 12;
    int next_z = z + ((i+1)>=12);
    if(next_z>=64)
        return true;
    if(fill[i][z] == -1){
        for(int t=0;t<64;t++){
            int node_comp = node_companion[i][t];
            int comp = z_companion[i][t];
            //printf("%d %d %d %d %d %d\n",i,z,t,node_comp,comp,pointer[node_comp]);
            if(!used[node_comp][comp] && !used[i][t] && pointer[node_comp]<(z+1+(z%2==0))){
                if(i!=node_comp) {
                    used[node_comp][comp] = used[i][t] = 1;
                    fill[i][z] = t;
                    fill[node_comp][pointer[node_comp]] = comp;
                    pointer[i]++;
                    pointer[node_comp]++;
                    if (dfs(next_i, next_z))
                        return true;
                    pointer[i]--;
                    pointer[node_comp]--;
                    fill[i][z] = -1;
                    fill[node_comp][pointer[node_comp]] = -1;
                    used[node_comp][comp] = used[i][t] = 0;
                }else{
                    used[i][t] = 1;
                    fill[i][z] = t;
                    pointer[i]++;
                    if (dfs(next_i, next_z))
                        return true;
                    pointer[i]--;
                    fill[i][z] = -1;
                    used[i][t] = 0;
                }
            }
        }
        return false;
    }else return dfs(next_i,next_z);

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



    for(int i=0;i<n;i++)
        for(int z=0;z<stripe;z++)
            fill[i][z] = -1,used[i][z] = 0;

    for(int i=0;i<n;i++)
        pointer[i] = 0;
/*
    if(dfs(0,0)){

        for(int i=0;i<n;i++) {
            for (int z = 0; z < stripe; z++)
                printf("%d ", fill[i][z]);
            printf("\n");
        }

        for(int i=0;i<n;i++)
            for(int z=0;z<stripe;z++)
                reorder[i][fill[i][z]] = z;

        for(int i=0;i<n;i++){
            for(int pos=0;pos<stripe;pos++){
                int z = fill[i][pos];
                int node_comp = node_companion[i][z];
                int com = z_companion[i][z];
                int com_pos = reorder[node_comp][com];
                printf("%d,%d:%d,%d\n",i,pos,node_comp,com_pos);
            }
            printf("\n");
        }

    }else
        printf("No solution found");
        */

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
    assert(n % (n - k) == 0);
    gf_init();
    init_companion(n, k);
    init_theta(n, k);
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


    //encode_t *((encode_t **)data_collected)[MAX_NODE];


    //for (i = 0; i < n; i++)
    //    ((encode_t **)data_collected)[i] = (encode_t *) data_collected[i];

    uint8_t inv_matrix[MAX_NODE][MAX_NODE];
    uint8_t final_matrix[MAX_NODE][MAX_NODE];
    uint8_t u_final_matrix[MAX_NODE][MAX_NODE];

    memset(inv_matrix, 0, sizeof(inv_matrix));
    n = q * t;

    int j;

    for (i = 0; i < q; i++)
        inv_matrix[i][i] = 1;

    uint8_t matrix[q][q];
    for (i = 0; i < q; i++)
        for (j = 0; j < q; j++) {
            int id = y_0 * q + j;
            if (j == x_0)
                matrix[i][j] = theta[i][id];
            else
                matrix[i][j] = gf_mul(u, theta[i][id]);
        }

    for (i = 0; i < q; i++) {
        uint8_t f = matrix[i][i];

        if (!matrix[i][i]) {
            for (j = i + 1; j < q; j++)
                if (matrix[j][i])
                    break;
            assert(j != q);

            for (k = 0; k < q; k++) {
                uint8_t tmp = matrix[i][k];
                matrix[i][k] = matrix[j][k];
                matrix[j][k] = tmp;

                tmp = inv_matrix[i][k];
                inv_matrix[i][k] = inv_matrix[j][k];
                inv_matrix[j][k] = tmp;

            }

        }

        for (j = 0; j < q; j++) {
            matrix[i][j] = gf_div(matrix[i][j], f);
            inv_matrix[i][j] = gf_div(inv_matrix[i][j], f);
        }

        for (j = 0; j < q; j++)
            if (i != j) {
                f = matrix[j][i];
                for (k = 0; k < q; k++) {
                    matrix[j][k] ^= gf_mul(matrix[i][k], f);
                    inv_matrix[j][k] ^= gf_mul(inv_matrix[i][k], f);
                }
            }
    }

    memset(final_matrix, 0, sizeof(final_matrix));
    for (i = 0; i < q; i++)
        for (j = 0; j < n; j++) {
            for (k = 0; k < q; k++) {
                final_matrix[i][j] ^= gf_mul(inv_matrix[i][k], theta[k][j]);
            }

        }

    for (i = 0; i < q; i++)
        for (j = 0; j < n; j++) {
            u_final_matrix[i][j] = gf_mul(final_matrix[i][j], u);
        }


    int z_num[MAX_STRIPE];
    int z_pos[MAX_STRIPE];
    memset(z_pos, -1, sizeof(z_pos));
    for (int z_id = 0; z_id < beta; z_id++) {
        z_num[z_id] = (z_id / _pow(q, t - y_0 - 1) * q + x_0) * _pow(q, t - y_0 - 1) + z_id % _pow(q, t - y_0 - 1);
        z_pos[z_num[z_id]] = z_id;
    }

    int z_comp_pos[MAX_NODE][MAX_STRIPE];

    for (int i = 0; i < n; i++)
        for (int z = 0; z < stripe_size; z++)
            z_comp_pos[i][z] = z_pos[z_companion[i][z]];

    encode_t kappa[MAX_NODE][REGION_BLOCKS + 3];
    memset(kappa, 0, sizeof(kappa));
    for (index = 0; index < block_size; index++) {

        for (int z_id = 0; z_id < beta; z_id++) {
            //

            if (q == 1) {

                for (j = 0; j < n; j++) {

                    int z = z_num[z_id];
                    int z_index = (block_size * z_id + index) * REGION_BLOCKS;

                    int companion = node_companion[j][z];

                    if (data_collected[j]) {
                        if (companion != j && data_collected[companion]) {
                            int new_z = z_pos[z_companion[j][z]];
                            int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;


                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                encode_t a_cur = xor_region(((encode_t **)data_collected)[j][z_index + w],
                                                            multiply_region(((encode_t **)data_collected)[companion][new_z_index + w], u));


                                for (int e = 0; e < 1; e++)
                                    kappa[e][w] = xor_region(kappa[e][w], multiply_region(a_cur,
                                                                                          final_matrix[e][j]));

                            }


                        } else {
                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                for (int e = 0; e < 1; e++) {
                                    kappa[e][w] = xor_region(kappa[e][w], multiply_region(((encode_t **)data_collected)[j][z_index + w],
                                                                                          final_matrix[e][j]));
                                }

                            }

                        }
                    } else if (j != companion && data_collected[companion]) {

                        int new_z = z_pos[z_companion[j][z]];
                        int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            for (int e = 0; e < 1; e++) {
                                kappa[e][w] = xor_region(kappa[e][w],
                                                         multiply_region(((encode_t **)data_collected)[companion][
                                                                                 new_z_index + w],
                                                                         u_final_matrix[e][j]));
                            }

                        }
                    }
                }
            } else if (q == 2) {

                for (j = 0; j < n; j++) {

                    int z = z_num[z_id];
                    int z_index = (block_size * z_id + index) * REGION_BLOCKS;

                    int companion = node_companion[j][z];

                    if (data_collected[j]) {
                        if (companion != j && data_collected[companion]) {
                            int new_z = z_pos[z_companion[j][z]];
                            int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;


                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                encode_t a_cur = xor_region(((encode_t **)data_collected)[j][z_index + w],
                                                            multiply_region(((encode_t **)data_collected)[companion][new_z_index + w], u));


                                for (int e = 0; e < 2; e++)
                                    kappa[e][w] = xor_region(kappa[e][w], multiply_region(a_cur,
                                                                                          final_matrix[e][j]));
                            }


                        } else {
                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                for (int e = 0; e < 2; e++) {
                                    kappa[e][w] = xor_region(kappa[e][w], multiply_region(((encode_t **)data_collected)[j][z_index + w],
                                                                                          final_matrix[e][j]));
                                }

                            }

                        }
                    } else if (j != companion && data_collected[companion]) {

                        int new_z = z_pos[z_companion[j][z]];
                        int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            for (int e = 0; e < 2; e++) {
                                kappa[e][w] = xor_region(kappa[e][w],
                                                         multiply_region(((encode_t **)data_collected)[companion][
                                                                                 new_z_index + w],
                                                                         u_final_matrix[e][j]));
                            }


                        }
                    }
                }
            } else if (q == 3) {

                for (j = 0; j < n; j++) {

                    int z = z_num[z_id];
                    int z_index = (block_size * z_id + index) * REGION_BLOCKS;

                    int companion = node_companion[j][z];

                    if (data_collected[j]) {
                        if (companion != j && data_collected[companion]) {
                            int new_z = z_pos[z_companion[j][z]];
                            int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;


                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                encode_t a_cur = xor_region(((encode_t **)data_collected)[j][z_index + w],
                                                            multiply_region(((encode_t **)data_collected)[companion][new_z_index + w], u));


                                for (int e = 0; e < 3; e++)
                                    kappa[e][w] = xor_region(kappa[e][w], multiply_region(a_cur,
                                                                                          final_matrix[e][j]));


                            }


                        } else {
                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                for (int e = 0; e < 3; e++) {
                                    kappa[e][w] = xor_region(kappa[e][w], multiply_region(((encode_t **)data_collected)[j][z_index + w],
                                                                                          final_matrix[e][j]));
                                }



                            }

                        }
                    } else if (j != companion && data_collected[companion]) {

                        int new_z = z_pos[z_companion[j][z]];
                        int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            for (int e = 0; e < 3; e++) {
                                kappa[e][w] = xor_region(kappa[e][w],
                                                         multiply_region(((encode_t **)data_collected)[companion][
                                                                                 new_z_index + w],
                                                                         u_final_matrix[e][j]));
                            }



                        }
                    }
                }
            } else if (q == 4) {

                for (j = 0; j < n; j++) {

                    int z = z_num[z_id];
                    int z_index = (block_size * z_id + index) * REGION_BLOCKS;

                    int companion = node_companion[j][z];

                    if (data_collected[j]) {
                        if (companion != j && data_collected[companion]) {
                            int new_z = z_pos[z_companion[j][z]];
                            int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;


                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                encode_t a_cur = xor_region(((encode_t **)data_collected)[j][z_index + w],
                                                            multiply_region(((encode_t **)data_collected)[companion][new_z_index + w], u));


                                for (int e = 0; e < 4; e++)
                                    kappa[e][w] = xor_region(kappa[e][w], multiply_region(a_cur,
                                                                                          final_matrix[e][j]));

                            }


                        } else {
                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                for (int e = 0; e < 4; e++) {
                                    kappa[e][w] = xor_region(kappa[e][w], multiply_region(((encode_t **)data_collected)[j][z_index + w],
                                                                                          final_matrix[e][j]));
                                }
                            }

                        }
                    } else if (j != companion && data_collected[companion]) {

                        int new_z = z_pos[z_companion[j][z]];
                        int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            for (int e = 0; e < 4; e++) {
                                kappa[e][w] = xor_region(kappa[e][w],
                                                         multiply_region(((encode_t **)data_collected)[companion][new_z_index + w],
                                                                         u_final_matrix[e][j]));
                            }

                        }
                    }
                }
            } else {
                assert(0);
            }


            for (j = 0; j < q; j++) {

                int z = (z_id / _pow(q, t - y_0 - 1) * q + j) * _pow(q, t - y_0 - 1) + z_id % _pow(q, t - y_0 - 1);
                int z_index = (z * block_size + index) * REGION_BLOCKS;
                //if(_mm256_extract_epi8(kappa[j][0],0) == 0)
                //    printf("%d %d %d\n",z_id,z,j);

                for (int w = 0; w < REGION_BLOCKS; w++) {
                    _mm256_stream_si256((encode_t *) memory + z_index + w, kappa[j][w]);
                    kappa[j][w] = _mm256_setzero_si256();
                }
            }

        }
    }
}


uint8_t a = 71;
uint8_t b = 201;

static void inline
sequential_decode(int index, int block_size, int *errors, int error_cnt,
                  int *sigmas, int sigma_max, int *err_id, bool *is_error, int *ok, encode_t *data_collected[MAX_NODE],
                  encode_t tmp[MAX_NODE][MAX_STRIPE * REGION_BLOCKS + 3], uint8_t final_matrix[MAX_NODE][MAX_NODE],
                  int q, int t) {
    int n = q * t;
    int k = n - q;

    int stripe_size = _pow(q, t);

    int s = 0;
    int j = 0;

    //The construction of B.
    while (s <= sigma_max) {
        for (int z = 0; z < stripe_size; z++)
            if (sigmas[z] == s) {
                //Use four 'if' to enable loop unrolling, which can improve the performance for about 30%.
                if (error_cnt == 1) {
                    for (j = 0; j < k; j++) {
                        int node_id = ok[j];

                        int z_index = (block_size * z + index) * REGION_BLOCKS;

                        int companion = node_companion[node_id][z];

                        if (companion != node_id) {
                            int new_z = z_companion[node_id][z];
                            int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                encode_t a_cur = xor_region(((encode_t **)data_collected)[node_id][z_index + w],
                                                            multiply_region(((encode_t **)data_collected)[companion][new_z_index + w], u));


                                for (int e = 0; e < 1; e++)
                                    tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                               multiply_region(a_cur,
                                                                                               final_matrix[e][node_id]));


                                _mm_prefetch(&((encode_t **)data_collected)[ok[j + 1]][(block_size * z + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);

                                _mm_prefetch(&((encode_t **)data_collected)[node_companion[ok[j + 1]][z]][
                                                     (block_size * z_companion[ok[j + 1]][z] + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);


                            }

                        } else if (companion == node_id) {
                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                for (int e = 0; e < 1; e++) {

                                    tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                               multiply_region(((encode_t **)data_collected)[node_id][
                                                                                                       (block_size * z +
                                                                                                        index) *
                                                                                                       REGION_BLOCKS +
                                                                                                       w],
                                                                                               final_matrix[e][node_id]));

                                }
                                _mm_prefetch(&((encode_t **)data_collected)[ok[j + 1]][(block_size * z + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);
                                _mm_prefetch(&((encode_t **)data_collected)[node_companion[ok[j + 1]][z]][
                                                     (block_size * z_companion[ok[j + 1]][z] + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);

                            }

                        }
                    }
                } else if (error_cnt == 2) {
                    for (j = 0; j < k; j++) {
                        int node_id = ok[j];

                        int z_index = (block_size * z + index) * REGION_BLOCKS;

                        int companion = node_companion[node_id][z];

                        if (companion != node_id) {
                            int new_z = z_companion[node_id][z];
                            int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                encode_t a_cur = xor_region(((encode_t **)data_collected)[node_id][z_index + w],
                                                            multiply_region(((encode_t **)data_collected)[companion][new_z_index + w], u));


                                for (int e = 0; e < 2; e++)
                                    tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                               multiply_region(a_cur,
                                                                                               final_matrix[e][node_id]));


                                _mm_prefetch(&((encode_t **)data_collected)[ok[j + 1]][(block_size * z + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);

                                _mm_prefetch(&((encode_t **)data_collected)[node_companion[ok[j + 1]][z]][
                                                     (block_size * z_companion[ok[j + 1]][z] + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);


                            }

                        } else if (companion == node_id) {
                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                for (int e = 0; e < 2; e++) {

                                    tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                               multiply_region(((encode_t **)data_collected)[node_id][
                                                                                                       (block_size * z +
                                                                                                        index) *
                                                                                                       REGION_BLOCKS +
                                                                                                       w],
                                                                                               final_matrix[e][node_id]));

                                }
                                _mm_prefetch(&((encode_t **)data_collected)[ok[j + 1]][(block_size * z + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);
                                _mm_prefetch(&((encode_t **)data_collected)[node_companion[ok[j + 1]][z]][
                                                     (block_size * z_companion[ok[j + 1]][z] + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);

                            }

                        }
                    }

                } else if (error_cnt == 3) {
                    for (j = 0; j < k; j++) {
                        int node_id = ok[j];

                        int z_index = (block_size * z + index) * REGION_BLOCKS;

                        int companion = node_companion[node_id][z];

                        if (companion != node_id) {
                            int new_z = z_companion[node_id][z];
                            int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                encode_t a_cur = xor_region(((encode_t **)data_collected)[node_id][z_index + w],
                                                            multiply_region(((encode_t **)data_collected)[companion][new_z_index + w], u));


                                for (int e = 0; e < 3; e++)
                                    tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                               multiply_region(a_cur,
                                                                                               final_matrix[e][node_id]));


                                _mm_prefetch(&((encode_t **)data_collected)[ok[j + 1]][(block_size * z + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);

                                _mm_prefetch(&((encode_t **)data_collected)[node_companion[ok[j + 1]][z]][
                                                     (block_size * z_companion[ok[j + 1]][z] + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);


                            }

                        } else if (companion == node_id) {
                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                for (int e = 0; e < 3; e++) {

                                    tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                               multiply_region(((encode_t **)data_collected)[node_id][
                                                                                                       (block_size * z +
                                                                                                        index) *
                                                                                                       REGION_BLOCKS +
                                                                                                       w],
                                                                                               final_matrix[e][node_id]));

                                }
                                _mm_prefetch(&((encode_t **)data_collected)[ok[j + 1]][(block_size * z + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);
                                _mm_prefetch(&((encode_t **)data_collected)[node_companion[ok[j + 1]][z]][
                                                     (block_size * z_companion[ok[j + 1]][z] + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);

                            }

                        }
                    }

                } else if (error_cnt == 4) {
                    for (j = 0; j < k; j++) {
                        int node_id = ok[j];

                        int z_index = (block_size * z + index) * REGION_BLOCKS;

                        int companion = node_companion[node_id][z];

                        if (companion != node_id) {
                            int new_z = z_companion[node_id][z];
                            int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                encode_t a_cur = xor_region(((encode_t **)data_collected)[node_id][z_index + w],
                                                            multiply_region(((encode_t **)data_collected)[companion][new_z_index + w], u));


                                for (int e = 0; e < 4; e++)
                                    tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                               multiply_region(a_cur,
                                                                                               final_matrix[e][node_id]));


                                _mm_prefetch(&((encode_t **)data_collected)[ok[j + 1]][(block_size * z + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);

                                _mm_prefetch(&((encode_t **)data_collected)[node_companion[ok[j + 1]][z]][
                                                     (block_size * z_companion[ok[j + 1]][z] + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);
                            }

                        } else if (companion == node_id) {
                            for (int w = 0; w < REGION_BLOCKS; w++) {
                                for (int e = 0; e < 4; e++) {

                                    tmp[e][z * REGION_BLOCKS + w] = xor_region(tmp[e][z * REGION_BLOCKS + w],
                                                                               multiply_region(((encode_t **)data_collected)[node_id][
                                                                                                       (block_size * z +
                                                                                                        index) *
                                                                                                       REGION_BLOCKS +
                                                                                                       w],
                                                                                               final_matrix[e][node_id]));

                                }
                                _mm_prefetch(&((encode_t **)data_collected)[ok[j + 1]][(block_size * z + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);
                                _mm_prefetch(&((encode_t **)data_collected)[node_companion[ok[j + 1]][z]][
                                                     (block_size * z_companion[ok[j + 1]][z] + index) * REGION_BLOCKS + w],
                                             _MM_HINT_NTA);

                            }

                        }
                    }

                } else {
                    assert(0);
                }

                for (j = 0; j < error_cnt; j++) {
                    int error = errors[j];

                    int companion = node_companion[error][z];
                    if (error != companion && !is_error[companion]) {
                        int new_z = z_companion[error][z];
                        int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

                        for (int w = 0; w < REGION_BLOCKS; w++) {
                            tmp[j][z * REGION_BLOCKS + w] = xor_region(tmp[j][z * REGION_BLOCKS + w],
                                                                       multiply_region(((encode_t **)data_collected)[companion][
                                                                                               new_z_index + w],
                                                                                       u));
                            if (j != error_cnt - 1)
                                _mm_prefetch(&((encode_t **)data_collected)[node_companion[errors[j + 1]][z]][
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
                            tmp[j][z_index] = tmp[comp_id][new_z_index] = _mm256_setzero_si256();


                            _mm256_stream_si256(&((encode_t **)data_collected)[error][(block_size * z + index) * REGION_BLOCKS + w], a_cur);
                            _mm256_stream_si256(&((encode_t **)data_collected)[companion][(block_size * new_z + index) * REGION_BLOCKS + w],
                                                a_companion);
                        }
                    } else if (companion == error || !is_error[companion]) {
                        for (int w = 0; w < REGION_BLOCKS; w++) {

                            _mm256_stream_si256(&((encode_t **)data_collected)[error][(block_size * z + index) * REGION_BLOCKS + w],
                                                tmp[j][z * REGION_BLOCKS + w]);
                            tmp[j][z*REGION_BLOCKS + w] = _mm256_setzero_si256();
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

    assert(error_cnt <= 4);

    int err_id[MAX_NODE];
    bool is_error[MAX_NODE];
    int ok[MAX_NODE + 1];
    //encode_t *data_ptr[MAX_NODE];

    int errors[error_cnt];

    int index = 0;
    for (int i = 0; i < n; i++)
        if (!data[i])
            errors[index++] = i;

    memset(is_error, 0, sizeof(bool) * n);

    for (int i = 0; i < error_cnt; i++) {
        is_error[errors[i]] = 1;
    }

    index = 0;

    for (int i = 0; i < n; i++)
        if (is_error[i]) {
            data[i] = memory_allocated[index++];
        }

    int sigmas[stripe_size];

    for (int z = 0; z < stripe_size; z++)
        sigmas[z] = compute_sigma(errors, error_cnt, q, t, z);

    int block_size = len / stripe_size / REGION_SIZE;

    //for (int i = 0; i < n; i++)
    //    data_ptr[i] = (encode_t *) data[i];

    encode_t tmp[MAX_NODE][MAX_STRIPE * REGION_BLOCKS + 3];
    for (int i = 0; i < error_cnt; i++)
        for (int z = 0; z < stripe_size * REGION_BLOCKS; z++)
            tmp[i][z] = _mm256_setzero_si256();

    for (int i = 0; i < error_cnt; i++)
        err_id[errors[i]] = i;


    ok[k] = 0;
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (!is_error[i])
            ok[j++] = i;
    }

    int sigma_max = 0;
    for (int z = 0; z < stripe_size; z++) {
        if (sigmas[z] > sigma_max)
            sigma_max = sigmas[z];
    }

    uint8_t inv_matrix[MAX_NODE][MAX_NODE];
    uint8_t final_matrix[MAX_NODE][MAX_NODE];
    uint8_t matrix[MAX_NODE][MAX_NODE];

    memset(inv_matrix, 0, sizeof(inv_matrix));

    for (int i = 0; i < error_cnt; i++)
        inv_matrix[i][i] = 1;


    memset(matrix, 0, sizeof(matrix));
    for (int i = 0; i < error_cnt; i++)
        for (j = 0; j < error_cnt; j++) {
            matrix[i][j] = theta[i][errors[j]];
        }


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

    for (index = 0; index < block_size; index++) {
        //memset(tmp,0,sizeof(tmp));
        sequential_decode(index, block_size, errors, error_cnt, sigmas, sigma_max, err_id, is_error, ok, ((encode_t **)data), tmp,
                          final_matrix, q, t);
    }

    //printf("Total Clock Time: %.2fs\n",(clock() - start)/(double)CLOCKS_PER_SEC);

    //printf("Encode Throughput: %.2fMB/s\n",(double)len * k/((clock() - start)/(double)CLOCKS_PER_SEC) * 1e-6 );
}

#define MAX_RS_NODE 64
/*
void rs_encode(int len,int row,int *errors, int k, int *alive, uint8_t **data,uint8_t **coding) {
    //data (k * len)
    //matrix (row * k)

    assert(len % sizeof(encode_t) == 0);
    int cur_row = 0;
    encode_t *((encode_t **)data_collected)[MAX_RS_NODE];
    encode_t *coding_ptr[MAX_RS_NODE];

    for (int i = 0; i < k; i++)
        ((encode_t **)data_collected)[i] = (encode_t *) data[i];

    for (int i = 0; i < row; i++)
        coding_ptr[i] = (encode_t *) coding[i];

    while (row >= 4) {
        for (int i = 0; i < len; i += sizeof(encode_t)) {
             for (int j = 0; j < k; j++)
                for (int r = 0; r < 4; r++) {
                    coding_ptr[cur_row + r][i] = xor_region(coding_ptr[cur_row + r][i],
                                                            multiply_region(((encode_t **)data_collected)[j][i], matrix[][]));
                }
        }
        row -= 4;
        cur_row += 4;
    }


}
 */
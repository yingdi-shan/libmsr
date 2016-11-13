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

typedef __m256i encode_t;

#define MAX_NODE 16
#define MAX_STRIPE 64
#define AVX2

#include <stdint.h>
#include <assert.h>
#include "gf.h"

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


int get_bit(int z, int y, int q, int t) {
    return z / _pow(q, y) % q;
}


int permute(int z, int remove_bit_id, int added_bit, int q, int t) {
    int result = 0;
    int power = _pow(q, t - 1);
    for (int i = t - 1; i >= 0; i--) {
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

#ifdef AVX2
__m256i mask_lo,mask_hi ;

__m256i low_table[256];
__m256i high_table[256];

void init_avx2(){
    mask_lo = _mm256_set1_epi8(0x0f);
    mask_hi = _mm256_set1_epi8(0xf0);
    //print_encode(mask_lo);
    //print_encode(mask_hi);

    int low_array[32];
    int high_array[32];
    for(int i=0;i<256;i++){
        for(int j=0;j<16;j++)
            low_array[j] = gf_mul(i,j);
        for(int j=0;j<16;j++)
            low_array[j + 16] = gf_mul(i,j);
        for(int j=0;j<16;j++)
            high_array[j] = gf_mul(i,j<<4);
        for(int j=0;j<16;j++)
            high_array[j + 16] = gf_mul(i,j<<4);
        low_table[i] = _mm256_setr_epi8(low_array[0],low_array[1],low_array[2],low_array[3],low_array[4],low_array[5],low_array[6],low_array[7],low_array[8],low_array[9],low_array[10],low_array[11],low_array[12],low_array[13],low_array[14],low_array[15],low_array[16],low_array[17],low_array[18],low_array[19],low_array[20],low_array[21],low_array[22],low_array[23],low_array[24],low_array[25],low_array[26],low_array[27],low_array[28],low_array[29],low_array[30],low_array[31]);
        high_table[i] = _mm256_setr_epi8(high_array[0],high_array[1],high_array[2],high_array[3],high_array[4],high_array[5],high_array[6],high_array[7],high_array[8],high_array[9],high_array[10],high_array[11],high_array[12],high_array[13],high_array[14],high_array[15],high_array[16],high_array[17],high_array[18],high_array[19],high_array[20],high_array[21],high_array[22],high_array[23],high_array[24],high_array[25],high_array[26],high_array[27],high_array[28],high_array[29],high_array[30],high_array[31]);

        //print_encode(low_table[i]);
        //print_encode(high_table[i]);
    }


}

inline __m256i xor_region(__m256i input1,__m256i input2){

    return _mm256_xor_si256(input1,input2);


}

inline __m256i multiply_region(__m256i input,uint8_t x){


    __m256i low = _mm256_and_si256(input,mask_lo);
    __m256i high = _mm256_and_si256(input,mask_hi);

    __m256i low_t = low_table[x],high_t = high_table[x];
    high = _mm256_srli_epi16(high,4);

    __m256i right = _mm256_shuffle_epi8(low_t,low);
    __m256i left = _mm256_shuffle_epi8(high_t,high);


    return _mm256_xor_si256(left,right);


}

void print_encode(__m256i x){
    //for(int i=0;i<32;i++)
    printf("%0x ",_mm256_extract_epi8(x,0));
    printf("%0x ",_mm256_extract_epi8(x,1));
    printf("%0x ",_mm256_extract_epi8(x,2));
    printf("%0x ",_mm256_extract_epi8(x,3));
    printf("%0x ",_mm256_extract_epi8(x,4));
    printf("%0x ",_mm256_extract_epi8(x,5));
    printf("%0x ",_mm256_extract_epi8(x,6));
    printf("%0x ",_mm256_extract_epi8(x,7));
    printf("%0x ",_mm256_extract_epi8(x,8));
    printf("%0x ",_mm256_extract_epi8(x,9));
    printf("%0x ",_mm256_extract_epi8(x,10));

    printf("\n");
}


#endif

static uint8_t node_companion[MAX_NODE][MAX_STRIPE];
static uint8_t z_companion[MAX_NODE][MAX_STRIPE];
static uint8_t theta[MAX_NODE][MAX_NODE];
static uint8_t u_theta[MAX_NODE][MAX_NODE];

static uint8_t inv_matrix[MAX_NODE][MAX_NODE];


//Used for symbol-remaping.
encode_t data_buffer[MAX_NODE][MAX_STRIPE];

//Input data chunk.
encode_t data_chunk[MAX_NODE][MAX_STRIPE];


static void invert_matrix(int *errors, int error_cnt) {

    memset(inv_matrix, 0, sizeof(inv_matrix));

    for (int i = 0; i < error_cnt; i++)
        inv_matrix[i][i] = 1;

    uint8_t matrix[error_cnt][error_cnt];
    for (int i = 0; i < error_cnt; i++)
        for (int j = 0; j < error_cnt; j++)
            matrix[i][j] = theta[i][errors[j]];

    for (int i = 0; i < error_cnt; i++) {
        uint8_t f = matrix[i][i];

        for (int j = 0; j < error_cnt; j++) {
            matrix[i][j] = gf_div(matrix[i][j], f);
            inv_matrix[i][j] = gf_div(inv_matrix[i][j], f);
        }
        //kappa[i] = gf_div(kappa[i], f);

        for (int j = 0; j < error_cnt; j++)
            if (i != j) {
                f = matrix[j][i];
                for (int k = 0; k < error_cnt; k++) {
                    matrix[j][k] ^= gf_mul(matrix[i][k], f);
                    inv_matrix[j][k] ^= gf_mul(inv_matrix[i][k], f);
                    //kappa[j] ^= gf_mul(kappa[i], f);
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


static void
compute_kappa(uint32_t stripe_size, int q, int t, int z, bool *errors, int error_cnt,
              uint8_t *minus_kappa) {

    for (int j = 0; j < error_cnt * sizeof(encode_t); j++)
        minus_kappa[j] = 0;

    int n = q * t;

    int k = n - q;


    for (int i = 0; i < k; i++)
        if (!errors[i]) {
            for (int j = 0; j < error_cnt; j++)
                for (int w = 0; w < sizeof(encode_t); w++) {

                    minus_kappa[j * sizeof(encode_t) + w] ^= gf_mul(theta[j][i], ((uint8_t *) data_chunk[i])[z * sizeof(encode_t) + w]);
                }
        }

    for (int i = k; i < k + error_cnt; i++)
        if (!errors[i])
            for (int w = 0; w < sizeof(encode_t); w++)
                minus_kappa[(i - k)*sizeof(encode_t) + w] ^= ((uint8_t *) data_chunk[i])[z * sizeof(encode_t) + w];


    for (int i = 0; i < k; i++) {
        int companion = node_companion[i][z];

        if (i != companion && (!errors[companion] || !errors[i])) {
            int z_comp = z_companion[i][z];
            for (int j = 0; j < error_cnt; j++)
                for (int w = 0; w < sizeof(encode_t); w++) {
                    minus_kappa[j * sizeof(encode_t) + w] ^= gf_mul(u_theta[j][i],
                                             ((uint8_t *) data_chunk[companion])[z_comp * sizeof(encode_t) + w]);
                }
        }
    }

    for (int i = k; i < k + error_cnt; i++) {
        int companion = node_companion[i][z];

        if (i != companion && (!errors[companion] || !errors[i]))
            for (int w = 0; w < sizeof(encode_t); w++) {
                minus_kappa[(i - k)*sizeof(encode_t) + w] ^= gf_mul(u,
                                             ((uint8_t *) data_chunk[companion])[z_companion[i][z] * sizeof(encode_t) +
                                                                                 w]);
            }
    }

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
    memset(res,0,sizeof(res));
    uint8_t * res_ptr = (uint8_t *)res;

    for (int j = 0; j < error_cnt; j++) {
        //res[j] = 0;
        for (int k = 0; k < error_cnt; k++)
            for (int w = 0; w < sizeof(encode_t); w++)
                res_ptr[j * sizeof(encode_t) + w] ^= gf_mul(inv_matrix[j][k], kappa[k * sizeof(encode_t) + w]);
    }

    for (int j = 0; j < error_cnt; j++)
        ((encode_t *)kappa)[j] = res[j];

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


static void
systematic_encode(uint32_t stripe_size, int *errors, int error_cnt, int *sigmas,
                  int q, int t) {

    int n = q * t;
    int k = n - q;

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

    for (int j = 0; j < k; j++)
        for (int z = 0; z < stripe_size; z++)
            data_buffer[j][z] = data_chunk[j][z];



    //The construction of B.
    for (int z = 0; z < stripe_size; z++)

        for (int j = 0; j < k; j++) {


            int companion = node_companion[j][z];

            if (companion < j && !is_error[companion]) {
                int new_z = z_companion[j][z];

                encode_t a_cur = xor_region(data_chunk[j][z],multiply_region(data_chunk[companion][new_z],u));

                encode_t a_companion = xor_region(multiply_region(data_chunk[j][z],u),data_chunk[companion][new_z]);


                data_buffer[j][z ] = a_cur;
                data_buffer[companion][new_z] = a_companion;

            }
        }


    for (int z = 0; z < stripe_size; z++)
        for (int j = 0; j < error_cnt; j++) {


                encode_t res = _mm256_set1_epi64x(0);

                for (int i = 0; i < k; i++) {
                    res = xor_region(res,multiply_region(data_buffer[i][z],theta[errors[j] - k][i]));
                }
                data_buffer[errors[j]][z] = res;



        }


    while (s <= sigma_max) {


        for (int z = 0; z < stripe_size; z++)
            if (sigmas[z] == s) {

                for (int j = 0; j < error_cnt; j++)
                     {

                        int error = errors[j];

                        int companion = node_companion[error][z];
                        int new_z = z_companion[error][z];


                        if (companion < error && is_error[companion]) {



                            encode_t a_cur = xor_region(multiply_region(data_buffer[error][z],a),multiply_region(data_buffer[companion][new_z],b));

                            encode_t a_companion = xor_region(multiply_region(data_buffer[error][z],b),multiply_region(data_buffer[companion][new_z],a));


                            data_buffer[error][z] = a_cur;
                            data_buffer[companion][new_z] = a_companion;

                        }
                    }

            }


        s++;
    }

    for (int i = 0; i < error_cnt; i++)
        for (int j = 0; j < stripe_size; j++)
            data_chunk[errors[i]][j] = data_buffer[errors[i]][j];


}


static void
sequential_decode(uint32_t stripe_size, int *errors, int error_cnt, int *sigmas,
                  int q, int t) {

    //clock_t start = clock();
    //printf("Sequential Decode:\n");
    int n = q * t;
    int z_total = stripe_size;

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


/*
    for (int i = 0; i < n; i++)
        for (int z = 0; z < pow(q, t); z++)
            test_companion(i, z, q, t);
*/

    //matrix L^(-1) = (1/(1-u^2),-u/(1-u^2);-u/(1-u^2),1/(1-u^2))
    uint8_t a = gf_div(1, gf_mul(1 ^ u, 1 ^ u));
    uint8_t b = gf_div(u, gf_mul(1 ^ u, 1 ^ u));

    assert(a ^ gf_mul(b, u) == 1);
    assert(b ^ gf_mul(a, u) == u);

    int s = 0;

    uint8_t *data_chunk_ptr[n];

    for (int i = 0; i < n; i++) {
        data_chunk_ptr[i] = (uint8_t *) data_chunk[i];
    }

    while (s <= sigma_max) {
        for (int z = 0; z < z_total; z++)
            if (sigmas[z] == s) {
                encode_t kappa[error_cnt];

                compute_kappa(stripe_size, q, t, z, is_error, error_cnt, (uint8_t *)kappa);
                //printf("Time for compute kappa:%lf\n",(clock() - start)/(double)(CLOCKS_PER_SEC));


                solve_equation(errors, error_cnt,(uint8_t *) kappa);
                //printf("Time for solve equation:%lf\n",(clock() - start)/(double)(CLOCKS_PER_SEC));
                for (int j = 0; j < error_cnt; j++)
                    data_chunk[errors[j]][z] = kappa[j];
            }


        for (int z = 0; z < stripe_size; z++)
            if (sigmas[z] == s) {

                for (int j = 0; j < error_cnt; j++)
                    for (int w = 0; w < sizeof(encode_t); w++) {

                        int error = errors[j];

                        int companion = node_companion[error][z];
                        int new_z = z_companion[error][z];




                        if (companion < error && is_error[companion]) {

                            int z_pos = z * sizeof(encode_t) + w;
                            int new_z_pos = new_z * sizeof(encode_t) + w;


                            uint8_t a_cur = gf_mul(data_chunk_ptr[error][z_pos], a) ^
                                            gf_mul(data_chunk_ptr[companion][new_z_pos], b);

                            uint8_t a_companion = gf_mul(data_chunk_ptr[error][z_pos], b) ^
                                                  gf_mul(data_chunk_ptr[companion][new_z_pos], a);


                            data_chunk_ptr[error][z_pos] = a_cur;
                            data_chunk_ptr[companion][new_z_pos] = a_companion;

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

    assert(len % (stripe_size * sizeof(encode_t)) == 0);

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

    invert_matrix(errors, error_cnt);

    //is_systematic = false;
    int sigmas[stripe_size];

    for (int z = 0; z < stripe_size; z++)
        sigmas[z] = compute_sigma(errors, error_cnt, q, t, z);

    int block_size = len / stripe_size / sizeof(encode_t);

    encode_t *data_ptr;

    for (index = 0; index < block_size; index++) {
        //Need to be optimized.
        //This is the performance bottom-neck.
        for (int i = 0; i < n; i++)
            if (!is_error[i]) {
                data_ptr = (encode_t *) data[i];
                //printf("%0x\n",*(uint8_t *)data_ptr);
                for (int j = 0; j < stripe_size; j++) {
                    data_chunk[i][j] = data_ptr[block_size * j + index];
                }
            }


        if (is_systematic)
            systematic_encode(stripe_size, errors, error_cnt, sigmas, q, t);
        else
            sequential_decode(stripe_size, errors, error_cnt, sigmas, q, t);


        for (int i = 0; i < n; i++)
            if (is_error[i]) {
                data_ptr = (encode_t *) data[i];

                for (int j = 0; j < stripe_size; j++)
                    data_ptr[block_size * j + index] = data_chunk[i][j];
            }

    }

    //free(errors);
    //free(data_buffer);
}


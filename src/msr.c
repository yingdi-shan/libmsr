//
// Created by syd on 17-6-12.
//
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include "msr.h"
#include "gf.h"
#include "arch.h"

//Magic number for CL-MSR.
const uint8_t u = 3;
const uint8_t a = 71;
const uint8_t b = 201;

static int ipow(int base, int exp) {
    int result = 1;
    while (exp) {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

static int get_bit(int z, int y, int q, int t) {
    return z / ipow(q, t - y - 1) % q;
}


static int permute(int z, int remove_bit_id, int added_bit, int q, int t) {
    int result = 0;
    int power = ipow(q, t - 1);
    int i;
    for (i = 0; i < t; i++) {
        if (i != remove_bit_id)
            result = result * q + (z / power) % q;
        else
            result = result * q + added_bit;
        power /= q;
    }
    return result;
}

static void init_companion(msr_conf *conf) {

    int i, z;

    for (i = 0; i < conf->n; i++)
        for (z = 0; z < conf->alpha; z++) {
            int x = i % conf->r;
            int y = i / conf->r;
            conf->node_companion[i * conf->alpha + z] = (uint8_t) (get_bit(z, y, conf->r, conf->groups) + y * conf->r);
            conf->z_companion[i * conf->alpha + z] = (uint8_t) (permute(z, y, x, conf->r, conf->groups));
        }
}

static void init_theta(msr_conf *conf) {
    int i, j;

    memset(conf->theta, 0, sizeof(uint8_t) * conf->r * conf->n);

    for (i = 0; i < conf->r; i++)
        for (j = 0; j < conf->k; j++) {
            conf->theta[i * conf->n + j] = gf_div(1, (uint8_t) (j ^ (i + conf->k)));
        }

    for (i = 0; i < conf->r; i++) {
        conf->theta[i * conf->n + i + conf->k] = 1;
    }
}

static void inline inverse_matrix(uint8_t *matrix, int n, uint8_t *inv) {
    int i=0, j=0;

    memset(inv, 0, n * n * sizeof(uint8_t));

    for (i = 0; i < n; i++)
        inv[i * n + i] = 1;

    for (i = 0; i < n; i++) {
        if (!matrix[i * n + i]) {
            for (j = i + 1; j < n; j++)
                if (matrix[j * n + i])
                    break;
            assert(j != n);

            for (int t = 0; t < n; t++) {
                uint8_t tmp = matrix[i * n + t];
                matrix[i * n + t] = matrix[j * n + t];
                matrix[j * n + t] = tmp;

                tmp = inv[i * n + t];
                inv[i * n + t] = inv[j * n + t];
                inv[j * n + t] = tmp;

            }
        }

        uint8_t f = matrix[i * n + i];

        for (j = 0; j < n; j++) {
            matrix[i * n + j] = gf_div(matrix[i * n + j], f);
            inv[i * n + j] = gf_div(inv[i * n + j], f);
        }

        for (j = 0; j < n; j++)
            if (i != j) {
                f = matrix[j * n + i];
                for (int t = 0; t < n; t++) {
                    matrix[j * n + t] ^= gf_mul(matrix[i * n + t], f);
                    inv[j * n + t] ^= gf_mul(inv[i * n + t], f);
                }
            }
    }
}

static int compute_sigma(int z, int *errors, int error_cnt, int r, int groups) {
    int sigma = 0;
    for (int i = 0; i < error_cnt; i++) {
        if (errors[i] % r == get_bit(z, errors[i] / r, r, groups)) {
            sigma++;
        }
    }
    return sigma;
}

int msr_init(msr_conf *conf, int n, int k,void* (*allocator)(size_t),void (*deallocator)(void *)) {
    int r = n - k;

    if (r <= 1 || r > k || n >= GF_SIZE - 1) {
        return -1;
    }

    gf_init();
    init_arch();

    conf->n = n;
    conf->k = k;
    conf->r = r;

    conf->groups = (n + r - 1) / r;
    conf->alpha = ipow(r, conf->groups);
    conf->beta = conf->alpha / r;

    conf->coding_unit_size = conf->alpha * REGION_SIZE * sizeof(uint8_t);

    conf->allocate = allocator;
    conf->deallocate = deallocator;

    conf->node_companion = allocator(sizeof(uint8_t) * conf->n * conf->alpha);
    conf->z_companion = allocator(sizeof(uint8_t) * conf->n * conf->alpha);
    conf->theta = allocator(sizeof(uint8_t) * conf->r * conf->n);

    init_companion(conf);
    init_theta(conf);

    return 0;
}

void msr_fill_encode_matrix(msr_encode_matrix *matrix, const msr_conf *conf, uint8_t **survived){
    int survived_cnt = 0;
    int erased_cnt = 0;
    int i,j,z;

    matrix->is_erased = conf->allocate(conf->n * sizeof(bool));
    matrix->erase_id = conf->allocate(conf->n * sizeof(int));

    for(i=0;i<conf->n;i++)
        if(survived[i] == NULL)
            matrix->is_erased[i] = true, erased_cnt++;
        else
            matrix->is_erased[i] = false, survived_cnt ++;

    matrix->survive_cnt = survived_cnt;
    matrix->erase_cnt = erased_cnt;

    matrix->survived = conf->allocate(survived_cnt * sizeof(int) + 1);
    matrix->erased = conf->allocate(erased_cnt * sizeof(int) + 1);

    survived_cnt = erased_cnt = 0;

    for(i=0;i<conf->n;i++)
        if(survived[i] == NULL) {
            matrix->erase_id[i] = erased_cnt;
            matrix->erased[erased_cnt++] = i;
        }
        else
            matrix->survived[survived_cnt++] = i;

    //Trick to prevent from overflow.
    matrix->survived[survived_cnt] = matrix->erased[erased_cnt] = 0;

    matrix->sigmas = conf->allocate(conf->alpha * sizeof(int));
    matrix->sigma_max = 0;
    for(z=0;z<conf->alpha;z++) {
        matrix->sigmas[z] = compute_sigma(z, matrix->erased, matrix->erase_cnt, conf->r, conf->groups);
        if (matrix->sigmas[z] > matrix->sigma_max)
            matrix->sigma_max = matrix->sigmas[z];
    }


    uint8_t * encode_matrix = conf->allocate(erased_cnt * erased_cnt * sizeof(uint8_t));
    uint8_t * inv_matrix = conf->allocate(erased_cnt * erased_cnt * sizeof(uint8_t));
    for (i = 0; i < erased_cnt; i++)
        for (j = 0; j < erased_cnt; j++) {
            encode_matrix[i * erased_cnt + j] = conf->theta[i * conf->n + matrix->erased[j]];
        }

    inverse_matrix(encode_matrix,erased_cnt,inv_matrix);

    matrix->matrix = conf->allocate(erased_cnt * conf->n * sizeof(uint8_t));

    memset(matrix->matrix, 0, erased_cnt * conf->n * sizeof(uint8_t));


    for (i = 0; i < erased_cnt; i++)
        for (j = 0; j < conf->n; j++) {
            for (int t = 0; t < erased_cnt; t++) {
                matrix->matrix[i * conf->n + j] ^= gf_mul(inv_matrix[i * erased_cnt + t], conf->theta[t * conf->n + j]);
            }
        }

    conf->deallocate(encode_matrix);
    conf->deallocate(inv_matrix);
}


static inline void decode_plane(const msr_encode_matrix* matrix, const msr_conf *conf, encode_t **data, encode_t *buf, size_t buf_len, int index, int block_size, int z) {
    int j;
    for(j = 0; j < conf->k;j++){
        int node_id = matrix->survived[j];
        int z_index = (block_size * z + index) * REGION_BLOCKS;
        int companion = conf->node_companion[node_id * conf->alpha + z];
        if (companion < conf->n && companion != node_id){
            int new_z = conf->z_companion[node_id * conf->alpha + z];
            int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;
            for(int w=0; w<REGION_BLOCKS;w++){
                encode_t a_cur = xor_region(data[node_id][z_index + w],multiply_region(data[companion][new_z_index + w],u));
                for(int e = 0; e < matrix->erase_cnt; e++){
                    buf[e * buf_len + z * REGION_BLOCKS + w] = xor_region(buf[e * buf_len + z * REGION_BLOCKS + w], multiply_region(a_cur,matrix->matrix[e * conf->n + node_id]));
                }
                prefetch(&data[matrix->survived[j + 1]][z_index + w]);
                prefetch(&data[conf->node_companion[matrix->survived[j+1] * conf->alpha + z]][(block_size * conf->z_companion[matrix->survived[j+1] + index]) * REGION_BLOCKS + w]);
            }
        } else {
            for(int w=0; w<REGION_BLOCKS;w++){
                for(int e = 0; e < matrix->erase_cnt; e++){
                    buf[e * buf_len + z * REGION_BLOCKS + w] = xor_region(buf[e * buf_len + z * REGION_BLOCKS + w],multiply_region(data[node_id][z_index + w],matrix->matrix[e * conf->n + node_id]));
                }
                prefetch(&data[matrix->survived[j + 1]][z_index + w]);
                prefetch(&data[conf->node_companion[matrix->survived[j+1] * conf->alpha + z]][(block_size * conf->z_companion[matrix->survived[j+1] + index]) * REGION_BLOCKS + w]);
            }
        }
    }

    for(j = 0; j < matrix->erase_cnt; j++) {
        int erased = matrix->erased[j];
        int companion = conf->node_companion[erased * conf->alpha + z];

        if(erased != companion && !matrix->is_erased[companion] && companion < conf->n){
            int new_z = conf->z_companion[erased * conf->alpha + z];
            int new_z_index = (block_size * new_z + index) * REGION_BLOCKS;

            for(int w = 0; w < REGION_BLOCKS; w++) {
                buf[j * buf_len + z * REGION_BLOCKS + w] = xor_region(buf[j * buf_len + z * REGION_BLOCKS + w],multiply_region(data[companion][new_z_index + w],u));
            }
        }
    }

}

static inline void write_plane(const msr_encode_matrix* matrix, const msr_conf *conf, encode_t **data,  encode_t *buf, size_t buf_len, int index, int block_size, int z) {
    int j;
    for(j = 0; j < matrix->erase_cnt; j++) {
        int erased = matrix->erased[j];

        int companion = conf->node_companion[erased * conf->alpha + z];
        int new_z = conf->z_companion[erased * conf->alpha + z];

        int comp_id = matrix->erase_id[companion];

        if(companion < erased && matrix->is_erased[companion]){
            for (int w = 0; w < REGION_BLOCKS; w++) {
                int z_index = z * REGION_BLOCKS + w;
                int new_z_index = new_z * REGION_BLOCKS + w;

                encode_t a_cur = xor_region(multiply_region(buf[j * buf_len + z_index],a),multiply_region(buf[comp_id * buf_len + new_z_index],b));
                encode_t a_companion = xor_region(multiply_region(buf[j * buf_len + z_index],b),multiply_region(buf[comp_id * buf_len + new_z_index],a));

                buf[j * buf_len + z_index] = buf[comp_id * buf_len + new_z_index] = zero();

                store(&data[erased][(block_size * z + index) * REGION_BLOCKS + w], a_cur);
                store(&data[companion][(block_size * new_z + index) * REGION_BLOCKS + w], a_companion);

            }
        } else if(companion == erased){
            for (int w = 0; w < REGION_BLOCKS; w++) {
                store(&data[erased][(block_size * z + index) * REGION_BLOCKS + w], buf[j * buf_len + z * REGION_BLOCKS + w]);
                buf[j * buf_len + z * REGION_BLOCKS + w] = zero();
            }
        }

    }
}

void msr_encode(int len, const msr_encode_matrix* matrix, const msr_conf *conf, uint8_t *buf, uint8_t **data, uint8_t **output){
    int z = 0,s = 0;

    int block_size = len / conf->alpha / REGION_SIZE;

    assert(len % (conf->alpha * REGION_SIZE) == 0);

    encode_t *input_ptr[conf->n];
    encode_t *buf_ptr = (encode_t *)buf;

    int erase_cnt = 0;
    for(int i=0; i < conf->n; i++){
        if(matrix->is_erased[i]){
            data[i] = output[erase_cnt ++];
        }
        input_ptr[i] = (encode_t *) data[i];
    }

    for(int index = 0; index < block_size; index++) {
        while (s <= matrix->sigma_max) {
            for (z = 0; z < conf->alpha; z++) {
                if (matrix->sigmas[z] == s) {
                    decode_plane(matrix, conf, input_ptr, buf_ptr, conf->alpha * REGION_BLOCKS, index, block_size, z);
                }
            }

            for (z = 0; z < conf->alpha; z++) {
                if (matrix->sigmas[z] == s) {
                    write_plane(matrix, conf, input_ptr , buf_ptr, conf->alpha * REGION_BLOCKS, index, block_size, z);
                }
            }
            s++;
        }
    }
}
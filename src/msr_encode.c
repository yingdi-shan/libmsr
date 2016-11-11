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


#define MAX_NODE 16
#define MAX_STRIPE 64

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


uint8_t node_companion[MAX_NODE][MAX_STRIPE];
uint8_t z_companion[MAX_NODE][MAX_STRIPE];
uint8_t theta[MAX_NODE][MAX_NODE];
uint8_t u_theta[MAX_NODE][MAX_NODE];

uint8_t inv_matrix[MAX_NODE][MAX_NODE];


void invert_matrix(int *errors, int error_cnt) {

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


void init_companion(int n, int k) {
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

void init_theta(int n, int k) {
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


void
compute_kappa(uint32_t stripe_size, uint8_t data_chunk[][stripe_size], int q, int t, int z, bool *errors, int error_cnt,
              uint8_t *minus_kappa) {

    for (int j = 0; j < error_cnt; j++)
        minus_kappa[j] = 0;

    int n = q * t;

    int k = n - q;

    for (int i = 0; i < k; i++)
        for (int j = 0; j < error_cnt; j++) {
            if (!errors[i])
                minus_kappa[j] ^= gf_mul(theta[j][i], data_chunk[i][z]);
        }

    for (int i = k; i < k + error_cnt; i++)
        if (!errors[i])
            minus_kappa[i - k] ^= data_chunk[i][z];

    for (int i = 0; i < k; i++) {
        int companion = node_companion[i][z];

        if (i != companion && (!errors[companion] || !errors[i])) {
            int z_comp = z_companion[i][z];
            for (int j = 0; j < error_cnt; j++) {
                minus_kappa[j] ^= gf_mul(u_theta[j][i], data_chunk[companion][z_comp]);
            }
        }
    }

    for (int i = k;i<k + error_cnt;i++){
        int companion = node_companion[i][z];

        if (i != companion && (!errors[companion] || !errors[i])) {
            minus_kappa[i-k] ^= gf_mul(u,data_chunk[companion][z_companion[i][z]]);
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


void solve_equation(int *errors, int error_cnt, uint8_t *kappa) {


    uint8_t res[error_cnt];

    for (int j = 0; j < error_cnt; j++) {
        res[j] = 0;
        for (int k = 0; k < error_cnt; k++)
            res[j] ^= gf_mul(inv_matrix[j][k], kappa[k]);
    }

    for (int j = 0; j < error_cnt; j++)
        kappa[j] = res[j];

}

//the number of each node range from 0 to n-1.
int compute_sigma(int *errors, int error_cnt, int q, int t, int z) {
    int sigma = 0;
    for (int i = 0; i < error_cnt; i++) {
        if (errors[i] % q == get_bit(z, errors[i] / q, q, t)) {
            sigma++;
        }
    }
    return sigma;
}


void systematic_encode(uint32_t stripe_size, uint8_t data_chunk[][stripe_size], int q, int t) {

    int n = q * t;



}


int
sequential_decode(uint32_t stripe_size, uint8_t data_chunk[][stripe_size], int *errors, int error_cnt, int q, int t) {

    //clock_t start = clock();

    int n = q * t;
    int z_total = stripe_size;
    int sigmas[z_total];
    int sigma_max = 0;
    for (int z = 0; z < z_total; z++) {
        sigmas[z] = compute_sigma(errors, error_cnt, q, t, z);

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


    while (s <= sigma_max) {
        for (int z = 0; z < z_total; z++)
            if (sigmas[z] == s) {
                uint8_t kappa[error_cnt];

                compute_kappa(stripe_size, data_chunk, q, t, z, is_error, error_cnt, kappa);
                //printf("Time for compute kappa:%lf\n",(clock() - start)/(double)(CLOCKS_PER_SEC));


                solve_equation(errors, error_cnt, kappa);
                //printf("Time for solve equation:%lf\n",(clock() - start)/(double)(CLOCKS_PER_SEC));
                for (int j = 0; j < error_cnt; j++)
                    data_chunk[errors[j]][z] = kappa[j];
            }

        char transformed[n][z_total];
        memset(transformed, 0, sizeof(char) * n * z_total);


        for (int z = 0; z < z_total; z++)
            if (sigmas[z] == s) {

                for (int j = 0; j < error_cnt; j++)
                    if (!transformed[errors[j]][z]) {

                        int error = errors[j];


                        int companion = node_companion[error][z];
                        int new_z = z_companion[error][z];


                        if (companion != error && is_error[companion]) {


                            //printf("%d %d\n",output_chunk[j][z],output_chunk[companion][permute(z, y, x, q, t)]);

                            uint8_t a_cur = gf_mul(data_chunk[error][z], a) ^
                                            gf_mul(data_chunk[companion][new_z], b);

                            uint8_t a_companion = gf_mul(data_chunk[error][z], b) ^
                                                  gf_mul(data_chunk[companion][new_z], a);


                            data_chunk[error][z] = a_cur;
                            data_chunk[companion][new_z] = a_companion;
                            transformed[error][z] = transformed[companion][new_z] = 1;
                        }
                    }
                //printf("Time for compute back:%lf\n",(clock() - start)/(double)(CLOCKS_PER_SEC));
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

    assert(len % (stripe_size) == 0);

    int error_cnt = 0;

    for (int i = 0; i < n; i++)
        if (!data[i])
            error_cnt++;

    uint8_t data_buffer[n][stripe_size];

    int errors[error_cnt];

    int index = 0;
    for (int i = 0; i < n; i++)
        if (!data[i])
            errors[index++] = i;

    bool is_error[n];

    memset(is_error, 0, sizeof(bool) * n);

    for (int i = 0; i < error_cnt; i++)
        is_error[errors[i]] = 1;

    index = 0;

    for (int i = 0; i < n; i++)
        if (is_error[i]) {
            data[i] = memory_allocated[index++];
        }

    invert_matrix(errors, error_cnt);

    for (int block = 0; block < len / (stripe_size); block++) {

        for (int i = 0; i < n; i++)
            if (!is_error[i]) {
                for (int j = 0; j < stripe_size; j++)
                    data_buffer[i][j] = data[i][len / stripe_size * j + block];
            }

        sequential_decode(stripe_size, data_buffer, errors, error_cnt, q, t);

        for (int i = 0; i < n; i++)
            if (is_error[i])
                for (int j = 0; j < stripe_size; j++)
                    data[i][len / stripe_size * j + block] = data_buffer[i][j];

    }

    //free(errors);
    //free(data_buffer);
}
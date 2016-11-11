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

#define INF 0x3fffffff

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

//Vandermonde matrix here, can be replaced by a cauchy matrix.
uint8_t get_theta(int row, int column) {
    return gf_pow(column + 1, row);
}

uint8_t *compute_kappa(uint32_t stripe_size,uint8_t data_chunk[][stripe_size], int q, int t, int z, bool *errors, int error_cnt) {
    uint8_t *minus_kappa = malloc(sizeof(uint8_t) * q);

    for (int j = 0; j < q; j++)
        minus_kappa[j] = 0;

    int n = q * t;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < error_cnt; j++) {

            if (!errors[i])
                minus_kappa[j] ^= gf_mul(get_theta(j, i), data_chunk[i][z]);
        }


    for (int i = 0; i < n; i++)
        for (int j = 0; j < error_cnt; j++) {
            int x = i % q;
            int y = i / q;

            //pos = 0;
            int z_y = get_bit(z, y, q, t);
            int companion = z_y + y * q;

            if (i != companion && !errors[companion]) {

                minus_kappa[j] ^= gf_mul(get_theta(j, i), gf_mul(u, data_chunk[companion][permute(z, y, x, q, t)]));
            } else if (i != companion && !errors[i]) {

                minus_kappa[j] ^= gf_mul(get_theta(j, i), gf_mul(u, data_chunk[companion][permute(z, y, x, q, t)]));
            }
        }




    //In GF(2^w), x = -x.
    return minus_kappa;
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
    uint8_t matrix[error_cnt][error_cnt];
    for (int i = 0; i < error_cnt; i++)
        for (int j = 0; j < error_cnt; j++)
            matrix[i][j] = get_theta(i, errors[j]);

    for (int i = 0; i < error_cnt; i++) {
        uint8_t f = matrix[i][i];
        for (int j = 0; j < error_cnt; j++) {
            matrix[i][j] = gf_div(matrix[i][j], f);
        }
        kappa[i] = gf_div(kappa[i], f);

        for (int j = 0; j < error_cnt; j++)
            if (i != j) {
                f = matrix[j][i];
                for (int k = 0; k < error_cnt; k++)
                    matrix[j][k] ^= gf_mul(matrix[i][k], f);
                kappa[j] ^= gf_mul(kappa[i], f);
            }
    }
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

/*
int in_set(int x, int *set, int set_cnt) {
    for (int i = 0; i < set_cnt; i++)
        if (set[i] == x)
            return 1;
    return 0;
}

int find_error_id(int x, int *errors, int error_cnt) {
    for (int i = 0; i < error_cnt; i++)
        if (errors[i] == x)
            return i;
    assert(0);
}
*/


/*  Sigma checker.
          int cnt = 0;
          int computed[12] = {0};
          memset(computed,0,sizeof(computed));
          for(int j=0;j<error_cnt;j++){
              int x = errors[j] % q;
              int y = errors[j] / q;
              int z_y = get_bit(z, y, q, t);
              int companion = y * q + z_y;
              if(computed[companion])
                  continue;
              if(in_set(companion,errors,error_cnt))
                  cnt++;
              computed[companion] = 1;

          }
          assert(cnt == sigmas[z]);
          */


int sequential_decode(uint32_t stripe_size,uint8_t data_chunk[][stripe_size], int *errors, int error_cnt, int q, int t) {
    int n = q * t;
    int z_total = stripe_size;
    int sigmas[z_total];
    int sigma_max = 0;
    for (int z = 0; z < z_total; z++) {
        sigmas[z] = compute_sigma(errors, error_cnt, q, t, z);

        if (sigmas[z] > sigma_max)
            sigma_max = sigmas[z];
    }

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
                uint8_t *kappa = compute_kappa(stripe_size,data_chunk, q, t, z, is_error, error_cnt);
                solve_equation(errors, error_cnt, kappa);
                for (int j = 0; j < error_cnt; j++)
                    data_chunk[errors[j]][z] = kappa[j];
                free(kappa);
            }

        char transformed[n][z_total];
        memset(transformed, 0, sizeof(char) * n * z_total);

        for (int z = 0; z < z_total; z++)
            if (sigmas[z] == s) {
                for (int j = 0; j < error_cnt; j++)
                    if (!transformed[errors[j]][z]) {

                        int error = errors[j];

                        int x = error % q;
                        int y = error / q;
                        int z_y = get_bit(z, y, q, t);
                        int companion = y * q + z_y;
                        int companion_z = permute(z, y, x, q, t);


                        if (z_y != x && is_error[companion]) {

                            assert(s != 0);


                            assert(sigmas[permute(z, y, x, q, t)] == s);
                            assert(!transformed[companion][companion_z]);

                            //printf("%d %d\n",output_chunk[j][z],output_chunk[companion][permute(z, y, x, q, t)]);

                            uint8_t a_cur = gf_mul(data_chunk[error][z], a) ^
                                            gf_mul(data_chunk[companion][companion_z], b);

                            uint8_t a_companion = gf_mul(data_chunk[error][z], b) ^
                                                  gf_mul(data_chunk[companion][companion_z], a);


                            data_chunk[error][z] = a_cur;
                            data_chunk[companion][companion_z] = a_companion;
                            transformed[error][z] = transformed[companion][companion_z] = 1;
                        }
                    }
            }

        s++;
    }

}

void msr_encode(int len, int n, int k, uint8_t **data,uint8_t **memory_allocated) {
    int r = n - k;

    assert(n % r == 0);

    int q = r;
    int t = n / q;


    uint32_t stripe_size = _pow(q, t);

    assert(len % (stripe_size ) == 0);

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

    for (int block = 0; block < len / (stripe_size); block++) {

        for (int i = 0; i < n; i++)
            if (!is_error[i]) {
                for (int j = 0; j < stripe_size; j++)
                    data_buffer[i][j] = data[i][len / stripe_size * j + block];
            }

        sequential_decode(stripe_size,data_buffer, errors, error_cnt, q, t);

        for (int i = 0; i < n; i++)
            if (is_error[i])
                for (int j = 0; j < stripe_size; j++)
                    data[i][len / stripe_size * j + block] = data_buffer[i][j];

    }

    //free(errors);
    //free(data_buffer);
}
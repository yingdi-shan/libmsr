//
// Created by syd on 16-11-11.
//
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "gf.h"
#include "msr.h"
#include <stdio.h>
#include <malloc.h>
#include <mm_malloc.h>

#define STRIPE_SIZE (_pow(q,t) * q * (t-1))
#define REGION_SIZE 512
#define DATA_SIZE (1<<30)

#define TEST_LOOP (10)

int main(){

    int n = 12;
    int r = 4;
    int k = n - r;

    int q = r;
    int t = n / q;

    init(n,k);

    uint8_t *data[n] ;
    uint8_t *memory_pre_allocated[k];

    for(int i=0;i<n;i++)
        data[i] = NULL;

    for(int i=0;i<n;i++) {
        posix_memalign((void *)&(data[i]),64,sizeof(uint8_t) * DATA_SIZE / k);
        //Warm the memory up.
        memset(data[i], 0xaa, sizeof(uint8_t) * DATA_SIZE / k);
    }

    for(int i=0;i<r;i++) {
        posix_memalign((void *)&(memory_pre_allocated[i]),64,sizeof(uint8_t) * DATA_SIZE / k);
        //Warm the memory up.
        memset(memory_pre_allocated[i],0x00,sizeof(uint8_t) * DATA_SIZE / k);
    }


    clock_t start = clock();

    for(int loop=0;loop<TEST_LOOP;loop++) {
        for(int i=0;i<r;i++)
           data[i + k] = NULL;
        msr_encode(DATA_SIZE / k, n, k, data, memory_pre_allocated);
    }

    printf("Total Clock Time: %.2fs\n",(clock() - start)/(double)CLOCKS_PER_SEC);

    printf("Encode Throughput: %.2fMB/s\n",TEST_LOOP * (double)DATA_SIZE/((clock() - start)/(double)CLOCKS_PER_SEC) * 1e-6);


    start = clock();

    for(int loop=0;loop<TEST_LOOP;loop++) {
        for(int i=0;i<r;i++)
            data[i] = NULL;
        msr_encode(DATA_SIZE / k, n, k, data, memory_pre_allocated);
    }

    printf("Total Clock Time: %.2fs\n",(clock() - start)/(double)CLOCKS_PER_SEC);

    printf("Decode Throughput: %.2fMB/s\n",TEST_LOOP * (double)DATA_SIZE/((clock() - start)/(double)CLOCKS_PER_SEC) * 1e-6 );

    int error = 1;

    int y_0 = error / q;
    int x_0 = error % q;
    uint8_t *input[n];
    for (int j = 0; j < n; j++){
        if(j!=error){
            posix_memalign((void *)(&input[j]),64,sizeof(uint8_t) * DATA_SIZE / k);
            int total = _pow(q,t);
            int len = _pow(q,t-1);
            for(int z_id=0;z_id<len;z_id++){
                int z = (z_id / _pow(q, t - y_0 - 1) * q + x_0) * _pow(q, t - y_0 - 1) + z_id % _pow(q, t - y_0 - 1);
                memcpy(input[j] + z_id*(DATA_SIZE/k/total),data[j] + z*(DATA_SIZE/k/total),DATA_SIZE/k/total);
            }

        }else
            input[j] = NULL;
    }

    uint8_t * memory;

    posix_memalign((void *)(&memory),64,sizeof(uint8_t) * DATA_SIZE / k);
    memset(memory,0,DATA_SIZE/k);

    start = clock();
    for(int i=0;i<TEST_LOOP;i++)
    msr_regenerate(DATA_SIZE/k/q,n,k,input,memory);

    printf("Total Clock Time: %.2fs\n",(clock() - start)/(double)CLOCKS_PER_SEC);

    printf("Regenerate Throughput: %.2fMB/s\n",TEST_LOOP * (double)DATA_SIZE/k/((clock() - start)/(double)CLOCKS_PER_SEC) * 1e-6 );

    free(memory);


    return 0;
}

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
#define DATA_SIZE (1<<30)

#define TEST_LOOP 5

int main(){
    gf_init();

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

    for(int i=0;i<k;i++) {
        posix_memalign((void *)&(data[i]),64,sizeof(uint8_t) * DATA_SIZE / k);
        //data[i] = malloc(sizeof(uint8_t) * DATA_SIZE / k);
        //Handle Page Fault.
        memset(data[i], 0xaa, sizeof(uint8_t) * DATA_SIZE / k);
    }

    for(int i=0;i<r;i++) {
        posix_memalign((void *)&(memory_pre_allocated[i]),64,sizeof(uint8_t) * DATA_SIZE / k);
        //Handle Page fault
        memset(memory_pre_allocated[i],0x00,sizeof(uint8_t) * DATA_SIZE / k);
    }


    clock_t start = clock();

    for(int loop=0;loop<TEST_LOOP;loop++) {
        for(int i=0;i<r;i++)
            data[i + k] = NULL;
        msr_encode(DATA_SIZE / k, n, k, data, memory_pre_allocated);
    }


    printf("Total Clock Time: %.2fs\n",(clock() - start)/(double)CLOCKS_PER_SEC);

    printf("Throughput: %.2fMB/s\n",TEST_LOOP * (double)DATA_SIZE/((clock() - start)/(double)CLOCKS_PER_SEC) * 1e-6 );


    return 0;
}
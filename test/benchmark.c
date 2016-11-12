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

#define STRIPE_SIZE (_pow(q,t) * q * (t-1))
#define DATA_SIZE ((1<<27))

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
        data[i] = malloc(sizeof(uint8_t) * DATA_SIZE / k);

        memset(data[i], 0xaa, sizeof(uint8_t) * DATA_SIZE / k);
    }

    for(int i=0;i<r;i++) {
        memory_pre_allocated[i] = malloc(sizeof(uint8_t) * DATA_SIZE / k);
        memset(memory_pre_allocated[i],0x00,sizeof(uint8_t) * DATA_SIZE / k);
    }

    clock_t start = clock();

    msr_encode(DATA_SIZE/k,n,k,data,memory_pre_allocated);

    printf("Throughput: %.2fMB/s\n",(double)DATA_SIZE/((clock() - start)/(double)CLOCKS_PER_SEC) * 1e-6);


    return 0;
}
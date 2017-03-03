//
// Created by syd on 16-10-23.
//
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "gf.h"
#include "msr.h"
#include <stdio.h>
#include <mm_malloc.h>


#define STRIPE_SIZE (_pow(q,t) * q * (t-1))
#define BLOCK_SIZE 512
#define DATA_SIZE (STRIPE_SIZE * (4096))

int main(int argc, char **argv) {
    srand(time(0));
    //srand(2);

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
        for(int j=0;j<DATA_SIZE/k;j++)
            data[i][j] = rand()%256;
       // memset(data[i], 0xaa, sizeof(uint8_t) * DATA_SIZE / k);
    }

    for(int i=0;i<r;i++) {
        posix_memalign((void *)&(memory_pre_allocated[i]),64,sizeof(uint8_t) * DATA_SIZE / k);
    }

    msr_encode(DATA_SIZE/k,n,k,data,memory_pre_allocated);

    for (int i = 0; i < n; i++) {
        printf("Original %d: ", i);
        for (int s = 0; s < 64; s++)
            printf("%x ", data[i][s]);
        printf("\n");
    }
    printf("-----------Begin to test decode-----------\n");
    int test_turn = 10;

    for (int i = 0; i < test_turn; i++) {
        printf("Turn %d:\n", i);

        uint8_t *input[n];
        for (int j = 0; j < n; j++)
            input[j] = NULL;

        int ok_cnt = 0;
        //int num[] = {1,2,3,4,5,6,7,8};
        while (ok_cnt < k) {
            //int ok_id = num[ok_cnt];
            int ok_id = rand() % n;
            if (!input[ok_id]) {
                input[ok_id] = data[ok_id];
                ok_cnt++;
            }
        }

        for(int i=0;i<r;i++) {
            posix_memalign((void *)&(memory_pre_allocated[i]),64,sizeof(uint8_t) * DATA_SIZE / k);
            memset(memory_pre_allocated[i],0x00,sizeof(uint8_t) * DATA_SIZE / k);
        }


        msr_encode(DATA_SIZE/k,n,k,input,memory_pre_allocated);


        for(int j=0;j<n;j++)
            assert(!memcmp(input[j],data[j],DATA_SIZE/k));

        printf("Check OK!\n");

        for(int i=0;i<r;i++) {
            free(memory_pre_allocated[i]);
        }


    }

    printf("-----------Begin to test regenerate-----------\n");
    test_turn = 12;
    for (int i = 0; i < test_turn; i++) {
        int id = i;
        int y_0 = id / q;
        int x_0 = id % q;
        uint8_t *input[n];
        for (int j = 0; j < n; j++){
            if(j!=i){
                posix_memalign((void *)(&input[j]),64,sizeof(uint8_t) * DATA_SIZE / k);
                int total = _pow(q,t);
                int len = _pow(q,t-1);
                for(int z_id=0;z_id<len;z_id++){
                    int z = (z_id / _pow(q, t - y_0 - 1) * q + x_0) * _pow(q, t - y_0 - 1) + z_id % _pow(q, t - y_0 - 1);
                    memcpy(input[j] + z_id*DATA_SIZE/k/total,data[j] + z*DATA_SIZE/k/total,DATA_SIZE/k/total);
                }

            }else
                input[j] = NULL;
        }

        uint8_t * memory;

        posix_memalign((void *)(&memory),64,sizeof(uint8_t) * DATA_SIZE / k);

        msr_regenerate(DATA_SIZE/k/q,n,k,input,memory);

        assert(!memcmp(memory,data[i],DATA_SIZE/k));

        printf("Check OK!\n");

        free(memory);
        for(int j=0;j<n;j++)
            if(j!=i)
                free(input[j]);

    }

    for(int i=0;i<n;i++)
        free(data[i]);

    return 0;

}
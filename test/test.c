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
#include "stdio.h"


#define q 4
#define t 3
#define STRIPE_SIZE (pow(q,t) * q * (t-1))
#define DATA_SIZE (STRIPE_SIZE * 1)

int main(int argc, char **argv) {
    //srand(time(0));
    srand(2);
    gf_init();

    int n = q * t;
    int k = q * (t - 1);

    int r = q;

    uint8_t *data = malloc(sizeof(uint8_t) * DATA_SIZE);
    memset(data, 0x01, DATA_SIZE);

    uint8_t **output = malloc(sizeof(uint8_t *) * q);
    for (int i = 0; i < q; i++)
        output[i] = malloc(DATA_SIZE / k);

    msr_encode(data, DATA_SIZE, output, n, k);


    for(int i=0;i<q;i++){
        printf("Encoded %d: ",i+k);
        for(int s=0;s<pow(q,t);s++)
            printf("%x ",output[i][s]);
        printf("\n");
    }

    int test_turn = 10;

    for (int i = 0; i < test_turn; i++) {
        printf("Turn:%d\n",i);

        uint8_t **input = malloc(sizeof(uint8_t *) * n);
        for(int j=0;j<n;j++)
            input[j] = NULL;

        int ok_cnt = 0;
        //int num[] = {0,3,4,6,7,8,10,11};
        while (ok_cnt < k) {
            //int ok_id = num[ok_cnt];
            int ok_id = rand()%n;
            if (!input[ok_id]) {
                if (ok_id < k) {
                    input[ok_id] = malloc(sizeof(uint8_t) * DATA_SIZE / k);
                    memset(input[ok_id],0x01,DATA_SIZE/k);
                } else
                    input[ok_id] = output[ok_id - k];

                ok_cnt++;
            }
        }

        uint8_t **output = malloc(sizeof(uint8_t *) *n);

        for(int j=0;j<r;j++)
            output[j] = malloc(sizeof(uint8_t) * DATA_SIZE / k);

        msr_decode(input,DATA_SIZE,output,n,k);

        int index=0;
        for(int j=0;j<n;j++)
            if(!input[j])
                input[j] = output[index ++];

        for(int j=0;j<n;j++)
        for(int r=0;r<q;r++)
            if(input[j] == output[r])
        {
            printf("%d after repaired: ",j);
            for(int s=0;s<pow(q,t);s++){
                printf("%0x ",output[r][s]);
            }
            printf("\n");
        }


    }


    return 0;

}
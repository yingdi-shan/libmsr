//
// Created by syd on 16-11-1.
//

#ifndef MSR_H
#define MSR_H

uint32_t _pow(uint32_t a, int b);

/**@brief Fill the unavailable data.
 * @param len Length of each block of data.
 * @param n The number of total blocks.
 * @param k The number of systematic blocks.
 * @param data Array of pointers to data buffer. If the pointer is NULL,the corresponding data will be generated.
 * @param memory_allocated Array of pointers to the pre-allocated memory.
 * @returns none
 */
void msr_encode(int len, int n, int k, uint8_t **data,uint8_t **memory_allocated);


//len should be the length of collected data.
int msr_regenerate(int len, int n, int k, uint8_t **data, uint8_t *output);

void init(int n,int k);

#endif //MSR_H

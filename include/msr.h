//
// Created by syd on 16-11-1.
//

#ifndef MSR_H
#define MSR_H

#include <inttypes.h>

/**@brief Fill the unavailable data.
 * @param len Length of each block of data. len should be at least 32KB for 8+4.
 * @param n The number of total blocks.
 * @param k The number of systematic blocks.
 * @param data Array of pointers to data buffer. If the pointer is NULL,the corresponding data will be generated.j
 * @param output Array of pointers to the output memory.
 * @returns none
 */
void msr_encode(int len, int n, int k, uint8_t **data,uint8_t **output);

/**@brief Regenerate the unavailable data.
 * @param len Length of each block of data. len should be the length of collected data.
 * @param n The number of total blocks.
 * @param k The number of systematic blocks.
 * @param data Array of pointers to data buffer. If the pointer is NULL,the corresponding data will be generated.
 * @param output The address to the regenerated data.
 * @returns none
 */
int msr_regenerate(int len, int n, int k, uint8_t **data, uint8_t *output);


int msr_regenerate_batch(int len,int n,int k,uint8_t *data, int broken, uint8_t *output);


void init(int n,int k);


#endif //MSR_H
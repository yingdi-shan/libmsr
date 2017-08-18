//
// Created by syd on 16-11-1.
//

#ifndef MSR_H
#define MSR_H

#include <inttypes.h>
#include <stdbool.h>

struct msr_conf_t{
    int n;
    int k;
    int r;

    int alpha;
    int beta;
    int groups;

    size_t coding_unit_size;

    void* (*allocate)(size_t);
    void (*deallocate)(void *);

    //Private, should not be modified.
    uint8_t *node_companion;
    uint8_t *z_companion;
    uint8_t *theta;
};

typedef struct msr_conf_t msr_conf;

struct msr_encode_matrix_t{

    uint8_t *matrix;

    int *survived;
    int survive_cnt;


    int *erased;
    int erase_cnt;

    bool *is_erased;
    int *erase_id;

    int *sigmas;
    int sigma_max;
};

typedef struct msr_encode_matrix_t msr_encode_matrix;


void msr_fill_encode_matrix(msr_encode_matrix *matrix, const msr_conf *conf, uint8_t **data);

/**@brief Fill the unavailable data.
 * @param len Length of each block of data. len should be at least 512 * alpha.
 * @param n The number of total blocks.
 * @param k The number of systematic blocks.
 * @param data Array of pointers to data buffer. If the pointer is NULL,the corresponding data will be generated.
 * @param output Array of pointers to the output memory.
 * @returns none
 */
void msr_encode(int len, const msr_encode_matrix* matrix, const msr_conf *conf, uint8_t *buf, uint8_t **data, uint8_t **output);

/**@brief Regenerate the unavailable data.
 * @param len Length of each block of data. len should be the length of collected data.
 * @param n The number of total blocks.
 * @param k The number of systematic blocks.
 * @param data Array of pointers to data buffer. If the pointer is NULL,the corresponding data will be generated.
 * @param output The address to the regenerated data.
 * @returns none
 */
void msr_regenerate(int len, msr_conf *conf, uint8_t **input, uint8_t *output);


int msr_init(msr_conf *conf,int n,int k,void* (*allocate)(size_t),void (*deallocate)(void *));


#endif //MSR_H
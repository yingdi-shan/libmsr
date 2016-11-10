//
// Created by syd on 16-11-1.
//

#ifndef LIBMBR_MSR_H
#define LIBMBR_MSR_H

int msr_encode(uint8_t *input,uint32_t data_size,uint8_t **output,int n,int k);
int msr_decode(uint8_t **data, uint32_t data_size, uint8_t **output_data,int n, int k);
int msr_regenerate(char *data);


#endif //LIBMBR_MSR_H

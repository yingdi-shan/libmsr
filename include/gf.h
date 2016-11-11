//
// Created by syd on 16-11-8.
//

#ifndef GF_H
#define GF_H

#include <inttypes.h>

//#define GF_BITS 8
#define GF_MOD 0x11D
#define GF_SIZE (256)

uint8_t gf_mul(uint32_t a,uint32_t b);
uint8_t gf_div(uint32_t a,uint32_t b);
void gf_init();

#endif //GF_H

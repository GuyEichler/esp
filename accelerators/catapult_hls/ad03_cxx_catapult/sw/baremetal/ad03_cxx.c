/* Copyright (c) 2011-2021 Columbia University, System Level Design Group */
/* SPDX-License-Identifier: Apache-2.0 */

#ifndef __riscv
#include <stdio.h>
#include <stdlib.h>
#endif

#include <fixed_point.h>
#include <math.h>

#include <esp_accelerator.h>
#include <esp_probe.h>

typedef int8_t token_t;

static unsigned DMA_WORD_PER_BEAT(unsigned _st)
{
    return (sizeof(void *) / _st);
}


#define SLD_AD03_CXX 0x053
#define DEV_NAME "sld,ad03_cxx_catapult"

/* <<--params-->> */
const int32_t batch = 1;
const int32_t mode = 0;

static unsigned in_words_adj;
static unsigned out_words_adj;
static unsigned in_len;
static unsigned out_len;
static unsigned in_size;
static unsigned out_size;
static unsigned out_offset;
static unsigned mem_size;

const int32_t size = 128;

/* Size of the contiguous chunks for scatter/gather */
#define CHUNK_SHIFT 20
#define CHUNK_SIZE BIT(CHUNK_SHIFT)
#define NCHUNK(_sz) ((_sz % CHUNK_SIZE == 0) ?		\
			(_sz / CHUNK_SIZE) :		\
			(_sz / CHUNK_SIZE) + 1)

/* User defined registers */
/* <<--regs-->> */
#define AD03_CXX_BATCH_REG 0x40
#define AD03_CXX_MODE_REG 0x44

static token_t raw_inputs[10][128] = {
{ -53, -2, -6, -4, -19, -17, -17, -22, -16, -19, -22, -27, -25, -25, -25, -31, -27, -37, -36, -35, -36, -34, -39, -37, -40, -40, -43, -43, -43, -42, -43, -38, -53, -4, -10, -1, -10, -17, -19, -18, -15, -26, -27, -24, -32, -25, -23, -36, -28, -35, -31, -37, -41, -42, -43, -35, -35, -37, -38, -39, -41, -41, -44, -35, -41, -2, -18, -13, -13, -9, -22, -21, -12, -24, -18, -23, -29, -26, -29, -37, -32, -28, -33, -34, -36, -38, -36, -40, -30, -37, -42, -43, -42, -38, -40, -38, -46, -2, -8, -15, -15, -13, -31, -29, -21, -27, -25, -27, -30, -29, -29, -32, -30, -39, -34, -34, -40, -40, -42, -38, -40, -41, -41, -42, -42, -42, -44, -38, },
{ -53, -4, -10, -1, -10, -17, -19, -18, -15, -26, -27, -24, -32, -25, -23, -36, -28, -35, -31, -37, -41, -42, -43, -35, -35, -37, -38, -39, -41, -41, -44, -35, -41, -2, -18, -13, -13, -9, -22, -21, -12, -24, -18, -23, -29, -26, -29, -37, -32, -28, -33, -34, -36, -38, -36, -40, -30, -37, -42, -43, -42, -38, -40, -38, -46, -2, -8, -15, -15, -13, -31, -29, -21, -27, -25, -27, -30, -29, -29, -32, -30, -39, -34, -34, -40, -40, -42, -38, -40, -41, -41, -42, -42, -42, -44, -38, -45, -3, -15, -5, -12, -16, -23, -22, -14, -25, -22, -27, -30, -28, -26, -32, -32, -38, -30, -39, -38, -36, -44, -39, -36, -36, -40, -44, -42, -43, -41, -37, },
{ -41, -2, -18, -13, -13, -9, -22, -21, -12, -24, -18, -23, -29, -26, -29, -37, -32, -28, -33, -34, -36, -38, -36, -40, -30, -37, -42, -43, -42, -38, -40, -38, -46, -2, -8, -15, -15, -13, -31, -29, -21, -27, -25, -27, -30, -29, -29, -32, -30, -39, -34, -34, -40, -40, -42, -38, -40, -41, -41, -42, -42, -42, -44, -38, -45, -3, -15, -5, -12, -16, -23, -22, -14, -25, -22, -27, -30, -28, -26, -32, -32, -38, -30, -39, -38, -36, -44, -39, -36, -36, -40, -44, -42, -43, -41, -37, -47, -8, -22, -8, -20, -9, -23, -25, -19, -23, -21, -23, -27, -25, -33, -32, -35, -34, -38, -35, -39, -40, -39, -39, -38, -38, -40, -40, -41, -44, -39, -38, },
{ -46, -2, -8, -15, -15, -13, -31, -29, -21, -27, -25, -27, -30, -29, -29, -32, -30, -39, -34, -34, -40, -40, -42, -38, -40, -41, -41, -42, -42, -42, -44, -38, -45, -3, -15, -5, -12, -16, -23, -22, -14, -25, -22, -27, -30, -28, -26, -32, -32, -38, -30, -39, -38, -36, -44, -39, -36, -36, -40, -44, -42, -43, -41, -37, -47, -8, -22, -8, -20, -9, -23, -25, -19, -23, -21, -23, -27, -25, -33, -32, -35, -34, -38, -35, -39, -40, -39, -39, -38, -38, -40, -40, -41, -44, -39, -38, -53, -8, -10, -17, -9, -18, -15, -23, -12, -21, -22, -34, -22, -23, -28, -30, -32, -37, -34, -36, -40, -40, -39, -38, -37, -39, -47, -45, -41, -41, -40, -38, },
{ -45, -3, -15, -5, -12, -16, -23, -22, -14, -25, -22, -27, -30, -28, -26, -32, -32, -38, -30, -39, -38, -36, -44, -39, -36, -36, -40, -44, -42, -43, -41, -37, -47, -8, -22, -8, -20, -9, -23, -25, -19, -23, -21, -23, -27, -25, -33, -32, -35, -34, -38, -35, -39, -40, -39, -39, -38, -38, -40, -40, -41, -44, -39, -38, -53, -8, -10, -17, -9, -18, -15, -23, -12, -21, -22, -34, -22, -23, -28, -30, -32, -37, -34, -36, -40, -40, -39, -38, -37, -39, -47, -45, -41, -41, -40, -38, -45, -3, -12, -9, -10, -17, -21, -20, -16, -22, -18, -21, -27, -29, -33, -33, -25, -32, -24, -30, -37, -35, -41, -40, -35, -37, -48, -38, -41, -42, -43, -40, },
{ -47, -8, -22, -8, -20, -9, -23, -25, -19, -23, -21, -23, -27, -25, -33, -32, -35, -34, -38, -35, -39, -40, -39, -39, -38, -38, -40, -40, -41, -44, -39, -38, -53, -8, -10, -17, -9, -18, -15, -23, -12, -21, -22, -34, -22, -23, -28, -30, -32, -37, -34, -36, -40, -40, -39, -38, -37, -39, -47, -45, -41, -41, -40, -38, -45, -3, -12, -9, -10, -17, -21, -20, -16, -22, -18, -21, -27, -29, -33, -33, -25, -32, -24, -30, -37, -35, -41, -40, -35, -37, -48, -38, -41, -42, -43, -40, -44, -5, -6, -5, -15, -15, -21, -17, -18, -21, -19, -33, -24, -28, -28, -33, -33, -35, -35, -40, -32, -40, -42, -38, -40, -39, -43, -37, -40, -42, -40, -38, },
{ -53, -8, -10, -17, -9, -18, -15, -23, -12, -21, -22, -34, -22, -23, -28, -30, -32, -37, -34, -36, -40, -40, -39, -38, -37, -39, -47, -45, -41, -41, -40, -38, -45, -3, -12, -9, -10, -17, -21, -20, -16, -22, -18, -21, -27, -29, -33, -33, -25, -32, -24, -30, -37, -35, -41, -40, -35, -37, -48, -38, -41, -42, -43, -40, -44, -5, -6, -5, -15, -15, -21, -17, -18, -21, -19, -33, -24, -28, -28, -33, -33, -35, -35, -40, -32, -40, -42, -38, -40, -39, -43, -37, -40, -42, -40, -38, -42, -7, -6, -13, -23, -15, -18, -15, -16, -21, -26, -24, -21, -26, -30, -32, -35, -37, -30, -31, -36, -39, -40, -38, -38, -37, -41, -38, -42, -41, -44, -39, },
{ -45, -3, -12, -9, -10, -17, -21, -20, -16, -22, -18, -21, -27, -29, -33, -33, -25, -32, -24, -30, -37, -35, -41, -40, -35, -37, -48, -38, -41, -42, -43, -40, -44, -5, -6, -5, -15, -15, -21, -17, -18, -21, -19, -33, -24, -28, -28, -33, -33, -35, -35, -40, -32, -40, -42, -38, -40, -39, -43, -37, -40, -42, -40, -38, -42, -7, -6, -13, -23, -15, -18, -15, -16, -21, -26, -24, -21, -26, -30, -32, -35, -37, -30, -31, -36, -39, -40, -38, -38, -37, -41, -38, -42, -41, -44, -39, -48, -17, -12, -5, -11, -7, -19, -25, -20, -11, -18, -25, -25, -27, -29, -36, -28, -31, -30, -32, -31, -38, -35, -42, -36, -41, -41, -40, -38, -43, -41, -40, },
{ -44, -5, -6, -5, -15, -15, -21, -17, -18, -21, -19, -33, -24, -28, -28, -33, -33, -35, -35, -40, -32, -40, -42, -38, -40, -39, -43, -37, -40, -42, -40, -38, -42, -7, -6, -13, -23, -15, -18, -15, -16, -21, -26, -24, -21, -26, -30, -32, -35, -37, -30, -31, -36, -39, -40, -38, -38, -37, -41, -38, -42, -41, -44, -39, -48, -17, -12, -5, -11, -7, -19, -25, -20, -11, -18, -25, -25, -27, -29, -36, -28, -31, -30, -32, -31, -38, -35, -42, -36, -41, -41, -40, -38, -43, -41, -40, -48, -7, -4, -4, -15, -15, -21, -16, -21, -10, -23, -26, -27, -23, -29, -33, -33, -32, -32, -35, -34, -40, -38, -40, -39, -39, -43, -38, -43, -44, -40, -38, },
{ -42, -7, -6, -13, -23, -15, -18, -15, -16, -21, -26, -24, -21, -26, -30, -32, -35, -37, -30, -31, -36, -39, -40, -38, -38, -37, -41, -38, -42, -41, -44, -39, -48, -17, -12, -5, -11, -7, -19, -25, -20, -11, -18, -25, -25, -27, -29, -36, -28, -31, -30, -32, -31, -38, -35, -42, -36, -41, -41, -40, -38, -43, -41, -40, -48, -7, -4, -4, -15, -15, -21, -16, -21, -10, -23, -26, -27, -23, -29, -33, -33, -32, -32, -35, -34, -40, -38, -40, -39, -39, -43, -38, -43, -44, -40, -38, -50, -7, -12, -8, -13, -17, -19, -24, -16, -18, -22, -27, -22, -27, -26, -26, -33, -35, -29, -34, -37, -37, -39, -41, -35, -38, -40, -37, -38, -43, -39, -40 }};

static token_t raw_outputs[10][128] = {
{ -44, -6, -7, -10, -14, -12, -18, -17, -15, -18, -19, -23, -23, -24, -25, -27, -27, -29, -28, -29, -30, -33, -33, -33, -33, -34, -34, -35, -37, -38, -38, -36, -44, -6, -8, -11, -15, -13, -17, -16, -16, -19, -19, -23, -23, -24, -25, -27, -28, -29, -29, -29, -32, -33, -33, -33, -33, -34, -35, -35, -37, -38, -38, -36, -44, -6, -7, -11, -14, -13, -17, -17, -16, -19, -20, -24, -23, -25, -27, -29, -28, -30, -30, -30, -31, -34, -34, -34, -34, -35, -35, -36, -38, -39, -38, -36, -45, -7, -8, -11, -14, -13, -17, -18, -16, -20, -19, -25, -24, -25, -27, -29, -29, -30, -30, -30, -32, -34, -34, -34, -34, -36, -36, -36, -38, -39, -39, -36 },
{ -42, -6, -8, -11, -14, -12, -18, -17, -16, -18, -18, -23, -24, -24, -25, -27, -27, -29, -28, -29, -30, -33, -33, -33, -33, -34, -34, -35, -37, -38, -38, -36, -42, -6, -9, -11, -15, -13, -17, -16, -16, -19, -19, -23, -24, -25, -26, -28, -28, -29, -29, -28, -32, -33, -33, -33, -33, -34, -35, -35, -37, -38, -39, -36, -42, -7, -8, -11, -15, -13, -17, -16, -16, -19, -20, -24, -24, -25, -27, -29, -29, -30, -29, -29, -31, -34, -34, -34, -34, -35, -35, -36, -39, -39, -38, -36, -43, -7, -9, -12, -15, -13, -18, -17, -16, -20, -19, -25, -25, -25, -27, -29, -29, -30, -30, -29, -31, -34, -34, -34, -34, -36, -36, -36, -38, -39, -39, -37 },
{ -40, -8, -10, -11, -15, -12, -17, -17, -16, -19, -19, -23, -24, -25, -25, -28, -28, -29, -28, -29, -31, -33, -33, -32, -32, -34, -34, -34, -36, -38, -37, -35, -40, -8, -11, -11, -15, -12, -17, -16, -16, -20, -19, -23, -24, -25, -25, -28, -28, -29, -29, -28, -32, -33, -33, -33, -32, -33, -34, -34, -36, -37, -37, -35, -40, -8, -10, -11, -15, -12, -17, -16, -16, -19, -19, -23, -23, -25, -26, -29, -28, -29, -29, -29, -31, -33, -33, -33, -33, -34, -34, -35, -37, -38, -37, -34, -40, -8, -11, -11, -15, -13, -17, -17, -16, -20, -19, -24, -24, -25, -25, -28, -28, -29, -29, -29, -31, -33, -33, -33, -33, -34, -35, -34, -36, -37, -37, -35 },
{ -43, -9, -11, -11, -16, -14, -19, -19, -17, -20, -20, -25, -26, -27, -27, -30, -29, -30, -30, -31, -32, -35, -36, -35, -34, -36, -36, -35, -37, -39, -39, -36, -43, -9, -11, -11, -16, -14, -18, -18, -17, -20, -20, -25, -25, -26, -27, -29, -29, -30, -30, -30, -33, -35, -35, -35, -34, -35, -36, -35, -37, -38, -39, -36, -43, -9, -10, -11, -15, -13, -17, -17, -16, -19, -20, -24, -24, -25, -26, -29, -28, -29, -29, -30, -31, -34, -34, -34, -34, -35, -36, -36, -38, -38, -38, -35, -43, -8, -10, -11, -15, -13, -17, -17, -16, -19, -19, -24, -24, -24, -25, -28, -27, -28, -28, -29, -30, -33, -34, -34, -33, -35, -36, -35, -37, -38, -38, -35 },
{ -44, -10, -11, -10, -15, -13, -17, -17, -16, -19, -18, -24, -24, -25, -26, -28, -28, -29, -29, -30, -31, -34, -35, -35, -34, -35, -36, -35, -37, -38, -38, -36, -44, -10, -12, -10, -15, -13, -16, -16, -15, -19, -18, -23, -23, -24, -25, -27, -28, -29, -29, -29, -32, -34, -34, -34, -33, -34, -36, -35, -37, -38, -38, -36, -44, -9, -11, -10, -15, -12, -15, -14, -15, -17, -18, -22, -22, -23, -24, -27, -26, -27, -28, -28, -30, -33, -33, -34, -33, -34, -35, -35, -37, -38, -37, -35, -44, -9, -10, -9, -14, -12, -15, -14, -14, -17, -17, -21, -21, -22, -23, -25, -25, -27, -27, -27, -29, -32, -33, -33, -32, -34, -35, -34, -36, -38, -38, -35 },
{ -44, -9, -11, -10, -15, -13, -17, -17, -16, -19, -18, -24, -24, -25, -26, -28, -28, -29, -29, -30, -31, -34, -35, -35, -34, -35, -36, -35, -37, -39, -38, -36, -44, -9, -11, -10, -15, -13, -16, -16, -16, -19, -18, -23, -23, -24, -25, -27, -27, -29, -29, -29, -32, -34, -34, -35, -33, -35, -36, -35, -37, -38, -38, -36, -44, -8, -10, -10, -14, -13, -15, -15, -15, -18, -18, -23, -22, -23, -24, -27, -27, -28, -28, -29, -30, -33, -34, -34, -33, -35, -36, -36, -38, -38, -38, -35, -44, -8, -10, -10, -14, -12, -15, -15, -14, -18, -17, -22, -22, -23, -24, -26, -26, -27, -28, -28, -30, -33, -33, -34, -32, -35, -35, -35, -37, -38, -38, -35 },
{ -43, -8, -8, -9, -14, -12, -17, -17, -15, -18, -18, -23, -23, -25, -25, -28, -28, -29, -29, -30, -31, -33, -35, -34, -34, -35, -35, -35, -37, -39, -38, -36, -44, -8, -8, -9, -14, -13, -16, -16, -14, -18, -18, -23, -23, -24, -25, -27, -28, -29, -29, -29, -32, -34, -34, -35, -33, -35, -36, -35, -38, -38, -38, -36, -44, -8, -8, -9, -14, -12, -15, -15, -14, -18, -19, -23, -22, -24, -25, -28, -27, -28, -29, -29, -31, -34, -34, -35, -34, -35, -36, -36, -38, -38, -38, -35, -44, -8, -8, -9, -13, -12, -15, -16, -14, -18, -18, -23, -22, -23, -24, -27, -27, -28, -28, -29, -30, -33, -34, -34, -33, -35, -36, -35, -37, -38, -38, -36 },
{ -41, -8, -9, -10, -14, -12, -17, -17, -16, -18, -18, -23, -23, -24, -24, -26, -27, -28, -28, -29, -30, -33, -34, -34, -33, -35, -35, -35, -37, -38, -38, -36, -41, -9, -9, -10, -14, -13, -16, -17, -16, -18, -18, -23, -23, -23, -24, -26, -27, -28, -28, -28, -32, -33, -34, -35, -33, -34, -36, -35, -37, -38, -38, -36, -42, -9, -8, -10, -14, -12, -16, -16, -15, -18, -18, -23, -22, -23, -24, -26, -26, -28, -28, -28, -30, -33, -34, -35, -33, -35, -36, -36, -38, -38, -37, -35, -42, -8, -9, -10, -13, -12, -16, -17, -15, -18, -18, -22, -22, -23, -24, -26, -26, -27, -27, -28, -30, -33, -34, -34, -33, -35, -36, -35, -37, -38, -38, -35 },
{ -42, -9, -8, -8, -13, -12, -16, -16, -14, -17, -18, -23, -23, -24, -24, -27, -27, -28, -28, -29, -31, -33, -34, -34, -33, -35, -35, -34, -36, -37, -37, -35, -42, -9, -8, -8, -13, -13, -15, -16, -14, -18, -17, -22, -22, -23, -24, -26, -27, -28, -28, -28, -31, -33, -33, -34, -33, -34, -35, -35, -37, -37, -37, -35, -42, -9, -8, -8, -13, -12, -14, -15, -13, -17, -17, -22, -21, -22, -23, -26, -26, -27, -28, -28, -30, -33, -33, -34, -33, -34, -35, -35, -37, -37, -37, -34, -42, -8, -8, -8, -12, -12, -14, -15, -13, -17, -17, -21, -21, -22, -22, -25, -25, -27, -27, -27, -29, -32, -33, -33, -32, -34, -35, -34, -36, -37, -37, -34 },
{ -44, -10, -9, -9, -13, -12, -16, -17, -15, -17, -17, -23, -22, -23, -24, -26, -26, -27, -28, -29, -30, -33, -34, -35, -33, -35, -35, -34, -36, -38, -37, -35, -44, -11, -9, -9, -13, -12, -16, -16, -15, -17, -17, -22, -22, -22, -23, -25, -26, -27, -28, -28, -31, -33, -34, -35, -32, -34, -36, -35, -36, -37, -38, -35, -44, -10, -8, -9, -13, -12, -15, -15, -15, -16, -17, -22, -21, -22, -23, -25, -25, -27, -28, -28, -30, -33, -33, -35, -33, -34, -35, -35, -37, -37, -37, -35, -44, -10, -9, -9, -13, -12, -15, -15, -14, -17, -17, -22, -21, -21, -23, -25, -25, -26, -27, -27, -29, -32, -33, -34, -32, -34, -35, -34, -36, -37, -37, -35 }};


unsigned abs(const unsigned input)
{
    return input < 0 ? -input : input;
}

unsigned allowed_error = 0;

static int validate_buf(token_t *out, token_t *gold)
{
	int i;
	int j;
	unsigned errors = 0;

#ifndef __riscv
	printf("  gold output data @%p\n", gold);
	printf("       output data @%p\n", out);
#else
	print_uart("  gold output data @"); print_uart_addr((uintptr_t) gold); print_uart("\n");
	print_uart("       output data @"); print_uart_addr((uintptr_t) out); print_uart("\n");
#endif

	for (i = 0; i < batch; i++) {
		for (j = 0; j < size; j++)
        {
            token_t gold_data = gold[i * out_words_adj + j];
            token_t out_data = out[i * out_words_adj + j];
            unsigned error_it = abs(gold_data - out_data);

			if (error_it > allowed_error)
            {
				errors++;
            }
#if 0
#ifndef __riscv
        	printf("  [%X, %X] data %lX / gold %lX\n", i, j, out_data_fxd, gold_data_fxd);
#else
        	print_uart("  [");
            print_uart_int((uint32_t)i);
            print_uart(",");
            print_uart_int((uint32_t)j);
            print_uart("] data ");
            print_uart_int64(out_data_fxd);
            print_uart(" / gold ");
            print_uart_int64(gold_data_fxd);
			if (error_it > allowed_error)
                print_uart(" *** ERROR ***");
            print_uart("\n");
#endif
#endif
        }
    }

#ifndef __riscv
	printf("  total errors %u\n", errors);
#else
	print_uart("  total errors "); print_uart_int(errors); print_uart("\n");
#endif

	return errors;
}

static void init_buf (token_t *inputs, token_t * gold_outputs)
{
	int i;
	int j;

#ifndef __riscv
	printf("  input data @%p\n", inputs);
#else
	print_uart("       input  data @"); print_uart_addr((uintptr_t) inputs); print_uart("\n");
#endif

	for (i = 0; i < batch; i++)
    {
		for (j = 0; j < size; j++)
        {
			inputs[i * in_words_adj + j] = raw_inputs[i][j];
        }
    }

#ifndef __riscv
	printf("  gold output data @%p\n", gold_outputs);
#else
	print_uart("  gold output data @"); print_uart_addr((uintptr_t) gold_outputs); print_uart("\n");
#endif

    for (i = 0; i < batch; i++) {
		for (j = 0; j < size; j++) {
			gold_outputs[i * out_words_adj + j] = raw_outputs[i][j];

#if 0
#ifndef __riscv
        	printf("  [%X, %X] init gold %lX\n", i, j, gold[i * out_words_adj + j]);
#else
        	print_uart("  [");
            print_uart_int((uint32_t)i);
            print_uart(",");
            print_uart_int((uint32_t)j);
            print_uart("] init gold ");
            print_uart_int64(gold[i * out_words_adj + j]);
            print_uart("\n");
#endif
#endif
        }
    }
}


int main(int argc, char * argv[])
{
	int i;
	int n;
	int ndev;
	struct esp_device *espdevs;
	struct esp_device *dev;
	unsigned done;
	unsigned **ptable;
	token_t *mem;
	token_t *gold;
	unsigned errors = 0;

	if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0) {
		in_words_adj = size;
		out_words_adj = size;
	} else {
		in_words_adj = round_up(size, DMA_WORD_PER_BEAT(sizeof(token_t)));
		out_words_adj = round_up(size, DMA_WORD_PER_BEAT(sizeof(token_t)));
	}

    in_len = in_words_adj * (batch);
	out_len = out_words_adj * (batch);
	in_size = in_len * sizeof(token_t);
	out_size = out_len * sizeof(token_t);
	out_offset  = in_len;
	mem_size = (out_offset * sizeof(token_t)) + out_size;

	// Search for the device
#ifndef __riscv
	printf("Scanning device tree... \n");
#else
	print_uart("Scanning device tree... \n");
#endif

	ndev = probe(&espdevs, VENDOR_SLD, SLD_AD03_CXX, DEV_NAME);

#ifndef __riscv
	printf("Found %d devices: %s\n", ndev, DEV_NAME);
#else
	print_uart("Found "); print_uart_int(ndev); print_uart(" devices: sld,ad03_cxx\n");
#endif

	if (ndev == 0) {
#ifndef __riscv
		printf("ad03_cxx not found\n");
#else
		print_uart("ad03_cxx not found\n");
#endif
		return 0;
	}


#if 1
	// Allocate memory
	gold = aligned_malloc(out_size);
	mem = aligned_malloc(mem_size);
#ifndef __riscv
	printf("  memory buffer base-address = %p\n", mem);
#else
	print_uart("  memory buffer base-address = "); print_uart_addr((uintptr_t) mem); print_uart("\n");
#endif
		// Alocate and populate page table
	ptable = aligned_malloc(NCHUNK(mem_size) * sizeof(unsigned *));
	for (i = 0; i < NCHUNK(mem_size); i++)
		ptable[i] = (unsigned *) &mem[i * (CHUNK_SIZE / sizeof(token_t))];
#ifndef __riscv
	printf("  ptable = %p\n", ptable);
	printf("  nchunk = %lu\n", NCHUNK(mem_size));
#else
	print_uart("  ptable = "); print_uart_addr((uintptr_t) ptable); print_uart("\n");
	print_uart("  nchunk = "); print_uart_int(NCHUNK(mem_size)); print_uart("\n");
#endif

#ifndef __riscv
	printf("  Generate input...\n");
#else
	print_uart("  Generate input...\n");
#endif

    init_buf(mem, gold);

#ifndef __riscv
	printf("  ... input ready!\n");
#else
	print_uart("  ... input ready!\n");
#endif

	// Pass common configuration parameters
#endif

	for (n = 0; n < ndev; n++) {

		dev = &espdevs[n];

		// Check DMA capabilities
		if (ioread32(dev, PT_NCHUNK_MAX_REG) == 0) {
#ifndef __riscv
			printf("  -> scatter-gather DMA is disabled. Abort.\n");
#else
			print_uart("  -> scatter-gather DMA is disabled. Abort.\n");
#endif
			return 0;
		}

		if (ioread32(dev, PT_NCHUNK_MAX_REG) < NCHUNK(mem_size)) {
#ifndef __riscv
			printf("  -> Not enough TLB entries available. Abort.\n");
#else
			print_uart("  -> Not enough TLB entries available. Abort.\n");
#endif
			return 0;
		}
#if 0
		// Allocate memory
		gold = aligned_malloc(out_size);
		mem = aligned_malloc(mem_size);
#ifndef __riscv
		printf("  memory buffer base-address = %p\n", mem);
#else
		print_uart("  memory buffer base-address = "); print_uart_addr((uintptr_t) mem); print_uart("\n");
#endif
		// Alocate and populate page table
		ptable = aligned_malloc(NCHUNK(mem_size) * sizeof(unsigned *));
		for (i = 0; i < NCHUNK(mem_size); i++)
			ptable[i] = (unsigned *) &mem[i * (CHUNK_SIZE / sizeof(token_t))];
#ifndef __riscv
		printf("  ptable = %p\n", ptable);
		printf("  nchunk = %lu\n", NCHUNK(mem_size));
#else
		print_uart("  ptable = "); print_uart_addr((uintptr_t) ptable); print_uart("\n");
		print_uart("  nchunk = "); print_uart_int(NCHUNK(mem_size)); print_uart("\n");
#endif

#ifndef __riscv
		printf("  Generate input...\n");
#else
		print_uart("  Generate input...\n");
#endif

        init_buf(mem, gold);

#ifndef __riscv
		printf("  ... input ready!\n");
#else
		print_uart("  ... input ready!\n");
#endif

		// Pass common configuration parameters
#endif
		iowrite32(dev, SELECT_REG, ioread32(dev, DEVID_REG));
		iowrite32(dev, COHERENCE_REG, ACC_COH_NONE);

#ifndef __sparc
		iowrite32(dev, PT_ADDRESS_REG, (unsigned long) ptable);
#else
		iowrite32(dev, PT_ADDRESS_REG, (unsigned) ptable);
#endif
		iowrite32(dev, PT_NCHUNK_REG, NCHUNK(mem_size));
		iowrite32(dev, PT_SHIFT_REG, CHUNK_SHIFT);

		// Use the following if input and output data are not allocated at the default offsets
		//iowrite32(dev, SRC_OFFSET_REG, 0x0);
		//iowrite32(dev, DST_OFFSET_REG, 0x0);

		// Pass accelerator-specific configuration parameters
		/* <<--regs-config-->> */
		iowrite32(dev, AD03_CXX_BATCH_REG, batch);
		iowrite32(dev, AD03_CXX_MODE_REG, mode);

		// Flush (customize coherence model here)
		esp_flush(ACC_COH_NONE);

		// Start accelerators
#ifndef __riscv
		printf("  Start...\n");
#else
		print_uart("  Start...\n");
#endif

		iowrite32(dev, CMD_REG, CMD_MASK_START);

        // Wait for completion
		done = 0;
		while (!done) {
			done = ioread32(dev, STATUS_REG);
			done &= STATUS_MASK_DONE;
		}
		iowrite32(dev, CMD_REG, 0x0);

#ifndef __riscv
		printf("  Done\n");
		printf("  validating...\n");
#else
		print_uart("  Done\n");
		print_uart("  validating...\n");
#endif

		/* Validation */
		errors = validate_buf(&mem[out_offset], gold);
#ifndef __riscv
		if (errors)
			printf("  ... FAIL\n");
		else
			printf("  ... PASS\n");
#else
		if (errors)
			print_uart("  ... FAIL\n");
		else
			print_uart("  ... PASS\n");
#endif

		aligned_free(ptable);
		aligned_free(mem);
		aligned_free(gold);
	}

#ifndef __riscv
	printf("DONE\n");
#else
    print_uart("DONE\n");
#endif

	return 0;
}

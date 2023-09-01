// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef INC_ESPACC_H
#define INC_ESPACC_H

#include "../inc/espacc_config.h"
#include <cstdio>

#include <ap_fixed.h>
#include <ap_int.h>
#include "hls_linear_algebra.h"

#define __round_mask(x, y) ((y)-1)
#define round_up(x, y) ((((x)-1) | __round_mask(x, y))+1)

// Data types and constants
#define VALUES_PER_WORD (DMA_SIZE / DATA_BITWIDTH)
#if ((SIZE_IN_CHUNK_DATA % VALUES_PER_WORD) == 0)
#define SIZE_IN_CHUNK (SIZE_IN_CHUNK_DATA / VALUES_PER_WORD)
#else
#define SIZE_IN_CHUNK (SIZE_IN_CHUNK_DATA / VALUES_PER_WORD + 1)
#endif
#if ((SIZE_OUT_CHUNK_DATA % VALUES_PER_WORD) == 0)
#define SIZE_OUT_CHUNK (SIZE_OUT_CHUNK_DATA / VALUES_PER_WORD)
#else
#define SIZE_OUT_CHUNK (SIZE_OUT_CHUNK_DATA / VALUES_PER_WORD + 1)
#endif

// data word
#if (IS_TYPE_FIXED_POINT ==  1)
typedef ap_fixed<DATA_BITWIDTH,DATA_BITWIDTH-FRAC_BITS> word_t;
#elif (IS_TYPE_UINT == 1)
typedef ap_uint<DATA_BITWIDTH> word_t;
#elif (IS_TYPE_FLOAT == 1)
#if (DATA_BITWIDTH == 32)
typedef float word_t;
#else
#error "Floating point word bitwidth not supported. Only 32 is supported."
#endif
#else // (IS_TYPE_INT == 1)
typedef ap_int<DATA_BITWIDTH> word_t;
#endif

typedef struct dma_word {
    word_t word[VALUES_PER_WORD];
} dma_word_t;

typedef word_t in_data_word;
typedef word_t out_data_word;

// Ctrl
typedef struct dma_info {
    ap_uint<32> index;
    ap_uint<32> length;
    ap_uint<32> size;
} dma_info_t;

// The 'size' variable of 'dma_info' indicates the bit-width of the words
// processed by the accelerator. Here are the encodings:
#define SIZE_BYTE   0
#define SIZE_HWORD  1
#define SIZE_WORD   2
#define SIZE_DWORD  3

#if (DATA_BITWIDTH == 8)
#define SIZE_WORD_T SIZE_BYTE
#elif (DATA_BITWIDTH == 16)
#define SIZE_WORD_T SIZE_HWORD
#elif (DATA_BITWIDTH == 32)
#define SIZE_WORD_T SIZE_WORD
#else // if (DATA_BITWIDTH == 64)
#define SIZE_WORD_T SIZE_DWORD
#endif

#define X_MAX 6
#define Z_MAX 164

void top(dma_word_t *out, dma_word_t *in1,
	/* <<--params-->> */
	 const unsigned conf_info_iter,
	 const unsigned conf_info_x_dim,
	 const unsigned conf_info_z_dim,
	 dma_info_t &load_ctrl, dma_info_t &store_ctrl);

void compute(word_t _inbuff[SIZE_IN_CHUNK_DATA],
	     word_t _outbuff[SIZE_OUT_CHUNK_DATA]);

template<unsigned DIM>
int inverse(word_t new_mat[DIM][DIM], word_t out[DIM][DIM], unsigned real_dim);

//Traits
struct TRAITS_A:
	hls::matrix_multiply_traits<hls::NoTranspose,hls::NoTranspose,
	X_MAX,X_MAX,X_MAX,X_MAX,word_t, word_t>
{
	static const int ARCH = 2;
	//static const int INNER_II = 9;
	//static const int UNROLL_FACTOR = 3;
};

struct TRAITS_B:
	hls::matrix_multiply_traits<hls::NoTranspose,hls::Transpose,
	X_MAX,X_MAX,X_MAX,X_MAX,word_t, word_t>
{
	static const int ARCH = 2;
	//static const int INNER_II = 9;
	//static const int UNROLL_FACTOR = 3;
};

struct TRAITS_C:
	hls::matrix_multiply_traits<hls::NoTranspose,hls::NoTranspose,
	Z_MAX,X_MAX,X_MAX,X_MAX,word_t, word_t>
{
	static const int ARCH = 0;
	//static const int INNER_II = 9;
	//static const int UNROLL_FACTOR = 3;
};

struct TRAITS_D:
	hls::matrix_multiply_traits<hls::NoTranspose,hls::Transpose,
	Z_MAX,X_MAX,Z_MAX,X_MAX,word_t, word_t>
{
	static const int ARCH = 0;
	//static const int INNER_II = 9;
	//static const int UNROLL_FACTOR = 3;
};

/* struct TRAITS_E: */
/* 	hls::matrix_multiply_traits<hls::NoTranspose,hls::NoTranspose, */
/* 	X_MAX,Z_MAX,Z_MAX,Z_MAX,word_t, word_t> */
/* { */
/* 	static const int ARCH = 0; */
/* 	//static const int INNER_II = 9; */
/* 	//static const int UNROLL_FACTOR = 3; */
/* }; */

struct TRAITS_E:
	hls::matrix_multiply_traits<hls::Transpose,hls::NoTranspose,
	Z_MAX,X_MAX,Z_MAX,Z_MAX,word_t, word_t>
{
	static const int ARCH = 0;
	//static const int INNER_II = 9;
	//static const int UNROLL_FACTOR = 3;
};

struct TRAITS_F:
	hls::matrix_multiply_traits<hls::NoTranspose,hls::Transpose,
	Z_MAX,X_MAX,1,X_MAX,word_t, word_t>
{
	static const int ARCH = 0;
	//static const int INNER_II = 9;
	//static const int UNROLL_FACTOR = 3;
};

struct TRAITS_G:
	hls::matrix_multiply_traits<hls::Transpose,hls::Transpose,
	Z_MAX,1,X_MAX,Z_MAX,word_t, word_t>
{
	static const int ARCH = 0;
	//static const int INNER_II = 9;
	//static const int UNROLL_FACTOR = 3;
};

struct TRAITS_H:
	hls::matrix_multiply_traits<hls::NoTranspose,hls::Transpose,
	1, X_MAX,X_MAX,X_MAX,word_t, word_t>
{
	static const int ARCH = 2;
	//static const int INNER_II = 9;
	//static const int UNROLL_FACTOR = 3;
};

struct TRAITS_I:
	hls::matrix_multiply_traits<hls::NoTranspose,hls::NoTranspose,
	X_MAX,Z_MAX,Z_MAX,X_MAX,word_t, word_t>
{
	static const int ARCH = 0;
	//static const int INNER_II = 9;
	//static const int UNROLL_FACTOR = 3;
};

/* struct TRAITS_J: */
/* 	hls::matrix_multiply_traits<hls::NoTranspose,hls::Transpose, */
/* 	X_MAX,X_MAX,Z_MAX,X_MAX,word_t, word_t> */
/* { */
/* 	static const int ARCH = 0; */
/* 	//static const int INNER_II = 9; */
/* 	//static const int UNROLL_FACTOR = 3; */
/* }; */

struct TRAITS_J:
	hls::matrix_multiply_traits<hls::NoTranspose,hls::Transpose,
	Z_MAX,X_MAX,X_MAX,X_MAX,word_t, word_t>
{
	static const int ARCH = 0;
	//static const int INNER_II = 9;
	//static const int UNROLL_FACTOR = 3;
};

//Inverse traits
typedef hls::qr_inverse_traits<Z_MAX, word_t, word_t> INV_CFG;

struct MY_CONFIG : INV_CFG {
	struct QRF_CONFIG :
		hls::qrf_traits<Z_MAX,
		Z_MAX,
		word_t,
		/* INV_CFG::InternalType> { */
		word_t> {
		static const int ARCH = 0;
		static const int UPDATE_II = 1;
	};
	struct BACK_SUB_CONFIG :
		hls::back_substitute_traits<Z_MAX,
		/* INV_CFG::InternalType, */
		/* INV_CFG::InternalType> { */
		word_t,
		word_t> {
                //static const int ARCH = 1;
		static const int INNER_II = 1;
		static const int DIAG_II = 1;
	};
	struct MULTIPLIER_CONFIG :
		hls::matrix_multiply_traits<hls::NoTranspose,
		/* hls::ConjugateTranspose, */
		hls::NoTranspose,
		Z_MAX,
		Z_MAX,
		Z_MAX,
		Z_MAX,
		/* INV_CFG::InternalType, */
		word_t,
		word_t> {
                static const int ARCH = 0;
                static const int INNER_II = 1;
                //static const int UNROLL_FACTOR = 4;
	};
};

#endif

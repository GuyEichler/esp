/* Copyright (c) 2011-2023 Columbia University, System Level Design Group */
/* SPDX-License-Identifier: Apache-2.0 */

#include <stdio.h>
#ifndef __riscv
#include <stdlib.h>
#endif

#include <esp_accelerator.h>
#include <esp_probe.h>
#include <fixed_point.h>

#include "A_array.h"
#include "H_array.h"
#include "W_array.h"
#include "Q_array.h"
#include "initial_state_array.h"
#include "measurements_array.h"
#include "prediction_array.h"
#include "real_array.h"
#include "P_array.h"

//typedef int32_t token_t;
typedef float token_t;

static unsigned DMA_WORD_PER_BEAT(unsigned _st)
{
        return (sizeof(void *) / _st);
}


#define SLD_KALMAN_FX 0x149
#define DEV_NAME "sld,kalman_fx_vivado"

#define STATES 6
#define NEURONS 164
#define TIME_STAMPS 10
#define CHUNKS 1
#define BATCHES TIME_STAMPS / CHUNKS

/* <<--params-->> */
const int32_t inv_reset = 0;
const int32_t inv_num = 1;
const int32_t chunks = CHUNKS;
const int32_t iter = BATCHES;
const int32_t x_dim = STATES;
const int32_t z_dim = NEURONS;

static unsigned in_words_adj;
static unsigned in_words_adj_z;
static unsigned out_words_adj;
static unsigned in_len;
static unsigned out_len;
static unsigned in_size;
static unsigned out_size;
static unsigned out_offset;
static unsigned mem_size;

/* Size of the contiguous chunks for scatter/gather */
#define CHUNK_SHIFT 20
#define CHUNK_SIZE BIT(CHUNK_SHIFT)
#define NCHUNK(_sz) ((_sz % CHUNK_SIZE == 0) ?		\
			(_sz / CHUNK_SIZE) :		\
			(_sz / CHUNK_SIZE) + 1)

/* User defined registers */
/* <<--regs-->> */
#define KALMAN_FX_INV_RESET_REG 0x54
#define KALMAN_FX_INV_NUM_REG 0x50
#define KALMAN_FX_CHUNKS_REG 0x4C
#define KALMAN_FX_ITER_REG 0x48
#define KALMAN_FX_X_DIM_REG 0x44
#define KALMAN_FX_Z_DIM_REG 0x40

/*Helper functions*/


/*End of helper functions*/

static int validate_buf(token_t *out, token_t *gold)
{
	int i;
	int j;
	unsigned errors = 0;
	float MSE = 0.0;

	for (i = 0; i < iter; i++)
		for (j = 0; j < (x_dim + x_dim * x_dim)*chunks; j++)
		{

			token_t gold_val = gold[i * out_words_adj + j];
			token_t acc_val = out[i * out_words_adj + j];

			int32_t acc_fixed = float_to_fixed32(acc_val, 16);
			int32_t gold_fixed = float_to_fixed32(gold_val, 16);

			token_t diff;
			if(gold_val < acc_val)
				diff = acc_val - gold_val;
			else
				diff = gold_val - acc_val;

			MSE += diff * diff;

			int32_t diff_fixed = float_to_fixed32(diff, 16);

			if(j%(x_dim + x_dim*x_dim) < x_dim)
				printf("NOT ERROR: X Accelerator value: %d Golden value: %d index: %u iter: %d diff: %d \n", acc_fixed, gold_fixed, i * out_words_adj + j, i, diff_fixed);
			/* else */
			/* 	printf("NOT ERROR: P Accelerator value: %d Golden value: %d index: %u iter: %d diff: %d \n", acc_fixed, gold_fixed, i * out_words_adj + j, i, diff_fixed); */

			if (gold[i * out_words_adj + j] != out[i * out_words_adj + j])
			{
				if(diff / gold_val > 0.1 || diff / acc_val > 0.1 || diff / gold_val < -0.1 || diff / acc_val < -0.1){
					/* if(j < x_dim) */
					/* 	printf("ERROR: X Accelerator value: %d Golden value: %d index: %u iter: %d diff: %d \n", acc_fixed, gold_fixed, i * out_words_adj + j, i, diff_fixed); */
					/* else */
					/* 	printf("ERROR: P Accelerator value: %d Golden value: %d index: %u iter: %d diff: %d \n", acc_fixed, gold_fixed, i * out_words_adj + j, i, diff_fixed); */

					errors++;
				}
			}
		}

	MSE /= ((x_dim + x_dim * x_dim) * TIME_STAMPS);
	int32_t MSE_fixed = float_to_fixed32(MSE, 16);
	printf("Output MSE: %d \n", MSE);

	printf("Sum of errors is %u \n", errors);

	return errors;
}


static void init_buf (token_t *in, token_t * gold)
{
	int i;
	int j;


	for(i = 0; i < iter; i++)
	{//z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim

		j = 0;

		if(i == 0) //only for first iteration
		{
			//X
			for(; j < x_dim; j++)
			{
				in[j] = (token_t) initial[j];
			}

			//P
			for(; j < x_dim + x_dim * x_dim; j++)
			{
				in[j] = (token_t) 0.0;
			}

			//F
			for(; j < x_dim + x_dim * x_dim * 2; j++)
			{
				in[j] = (token_t) A[j - (x_dim + x_dim * x_dim)];
			}

			//Q
			for(; j < x_dim + x_dim * x_dim * 3; j++)
			{
				in[j] = (token_t) W[j - (x_dim + x_dim * x_dim * 2)];
			}

			//R
			for(; j < x_dim + x_dim * x_dim * 3 + z_dim * z_dim; j++)
			{
				in[j] = (token_t) Q[j - (x_dim + x_dim * x_dim * 3)];
			}

			//H
			for(; j < x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim; j++)
			{
				in[j] = (token_t) H[j - (x_dim + x_dim * x_dim * 3 + z_dim * z_dim)];
				//printf("Value of H = %f \n", measurements[NEURONS * (i+1) + j]);
			}
		}

		unsigned base_index = (x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim);

		//Z
		if(i == 0)
			for(; j < base_index + z_dim * chunks; j++)
			{

				in[j] = (token_t) measurements[NEURONS * (i+1) + j - base_index];
				int32_t val = float_to_fixed32(measurements[NEURONS * (i+1) + j], 3);
				//printf("Value of Z = %d index %d \n", val, j);
				/* in[in_words_adj + (i-1) * in_words_adj_z + j] = (token_t) measurements[NEURONS * (i+1) + j]; */
			}
		else
			for(; j < z_dim * chunks; j++)
			{
				in[in_words_adj + i * in_words_adj_z + j] = (token_t) measurements[NEURONS * (i * chunks + 1) + j];
				int32_t val = float_to_fixed32(measurements[NEURONS * i + j], 3);
			// if(i == 3)
			//printf("Value of Z = %d index %d \n", val, i * in_words_adj + j);
			}

	}


	/* for (i = 0; i < iter; i++) */
	/* 	for (j = 0; j < x_dim + x_dim * x_dim; j++) */
	/* 		gold[i * out_words_adj + j] = (token_t) j; */
	for(i = 0; i < TIME_STAMPS; i++)
		for(j = 0; j < x_dim + x_dim * x_dim; j++)
		{
			if(j < x_dim)
				gold[i * out_words_adj/chunks + j] = (token_t) prediction[STATES * (i+1) + j];
			else
			{
				gold[i * out_words_adj/chunks + j] = (token_t) P_flat[i * x_dim * x_dim + (j - x_dim)];
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
	unsigned coherence;

	if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0) {
		in_words_adj = x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim;
		in_words_adj_z = z_dim * chunks;
		out_words_adj = (x_dim + x_dim * x_dim) * chunks;
	} else {
		in_words_adj = round_up(x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim, DMA_WORD_PER_BEAT(sizeof(token_t)));
		in_words_adj_z = round_up(z_dim * chunks, DMA_WORD_PER_BEAT(sizeof(token_t)));
		out_words_adj = round_up((x_dim + x_dim * x_dim) * chunks, DMA_WORD_PER_BEAT(sizeof(token_t)));
	}
	in_len = in_words_adj + in_words_adj_z * iter;
	out_len = out_words_adj * (iter);
	in_size = in_len * sizeof(token_t);
	out_size = out_len * sizeof(token_t);
	out_offset  = in_len;
	mem_size = (out_offset * sizeof(token_t)) + out_size;


	// Search for the device
	printf("Scanning device tree... \n");

	ndev = probe(&espdevs, VENDOR_SLD, SLD_KALMAN_FX, DEV_NAME);
	if (ndev == 0) {
		printf("kalman_fx not found\n");
		return 0;
	}

	for (n = 0; n < ndev; n++) {

		printf("**************** %s.%d ****************\n", DEV_NAME, n);

		dev = &espdevs[n];

		// Check DMA capabilities
		if (ioread32(dev, PT_NCHUNK_MAX_REG) == 0) {
			printf("  -> scatter-gather DMA is disabled. Abort.\n");
			return 0;
		}

		if (ioread32(dev, PT_NCHUNK_MAX_REG) < NCHUNK(mem_size)) {
			printf("  -> Not enough TLB entries available. Abort.\n");
			return 0;
		}

		// Allocate memory
		gold = aligned_malloc(out_size);
		mem = aligned_malloc(mem_size);
		printf("  memory buffer base-address = %p\n", mem);

		// Alocate and populate page table
		ptable = aligned_malloc(NCHUNK(mem_size) * sizeof(unsigned *));
		for (i = 0; i < NCHUNK(mem_size); i++)
			ptable[i] = (unsigned *) &mem[i * (CHUNK_SIZE / sizeof(token_t))];

		printf("  ptable = %p\n", ptable);
		printf("  nchunk = %lu\n", NCHUNK(mem_size));

#ifndef __riscv
		for (coherence = ACC_COH_NONE; coherence <= ACC_COH_RECALL; coherence++) {
#else
		{
			/* TODO: Restore full test once ESP caches are integrated */
			coherence = ACC_COH_NONE;
#endif
			printf("  --------------------\n");
			printf("  Generate input...\n");
			init_buf(mem, gold);

			// Pass common configuration parameters

			iowrite32(dev, SELECT_REG, ioread32(dev, DEVID_REG));
			iowrite32(dev, COHERENCE_REG, coherence);

#ifndef __sparc
			iowrite32(dev, PT_ADDRESS_REG, (unsigned long long) ptable);
#else
			iowrite32(dev, PT_ADDRESS_REG, (unsigned) ptable);
#endif
			iowrite32(dev, PT_NCHUNK_REG, NCHUNK(mem_size));
			iowrite32(dev, PT_SHIFT_REG, CHUNK_SHIFT);

			// Use the following if input and output data are not allocated at the default offsets
			iowrite32(dev, SRC_OFFSET_REG, 0x0);
			iowrite32(dev, DST_OFFSET_REG, 0x0);

			// Pass accelerator-specific configuration parameters
			/* <<--regs-config-->> */
			iowrite32(dev, KALMAN_FX_INV_RESET_REG, inv_reset);
			iowrite32(dev, KALMAN_FX_INV_NUM_REG, inv_num);
			iowrite32(dev, KALMAN_FX_CHUNKS_REG, chunks);
			iowrite32(dev, KALMAN_FX_ITER_REG, iter);
			iowrite32(dev, KALMAN_FX_X_DIM_REG, x_dim);
			iowrite32(dev, KALMAN_FX_Z_DIM_REG, z_dim);

			// Flush (customize coherence model here)
			esp_flush(coherence);

			// Start accelerators
			printf("  Start...\n");
			iowrite32(dev, CMD_REG, CMD_MASK_START);

			// Wait for completion
			done = 0;
			while (!done) {
				done = ioread32(dev, STATUS_REG);
				done &= STATUS_MASK_DONE;
			}
			iowrite32(dev, CMD_REG, 0x0);

			printf("  Done\n");
			printf("  validating...\n");

			/* Validation */
			errors = validate_buf(&mem[out_offset], gold);
			if (errors)
				printf("  ... FAIL\n");
			else
				printf("  ... PASS\n");
		}
		aligned_free(ptable);
		aligned_free(mem);
		aligned_free(gold);
	}

	return 0;
}

// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "libesp.h"
#include "cfg.h"

#include "A_array.h"
#include "H_array.h"
#include "W_array.h"
#include "Q_array.h"
#include "initial_state_array.h"
#include "measurements_array.h"
#include "prediction_array.h"
#include "real_array.h"
#include "P_array.h"

static unsigned in_words_adj;
static unsigned in_words_adj_z;
static unsigned out_words_adj;
static unsigned in_len;
static unsigned out_len;
static unsigned in_size;
static unsigned out_size;
static unsigned out_offset;
static unsigned size;

/* User-defined code */
static int validate_buffer(token_t *out, token_t *gold)
{
	int i;
	int j;
	unsigned errors = 0;
	float MSE = 0.0;

	for (i = 0; i < iter; i++)
		for (j = 0; j < x_dim + x_dim * x_dim; j++)
		{

			token_t gold_val = gold[i * out_words_adj + j];
			token_t acc_val = out[i * out_words_adj + j];

			token_t diff;
			if(gold_val < acc_val)
				diff = acc_val - gold_val;
			else
				diff = gold_val - acc_val;

			MSE += diff * diff;

			if(j < x_dim)
				printf("NO ERROR: X Accelerator value: %f Golden value: %f index: %d iter: %d diff: %f\n", acc_val, gold_val, i * out_words_adj + j, i, diff);
			/* else */
			/* 	printf("N ERROR: P Accelerator value: %f Golden value: %f index: %d\n iter: %d diff: %f", acc_val, gold_val, i * out_words_adj + j, i, diff); */

			if (gold[i * out_words_adj + j] != out[i * out_words_adj + j])
			{
				if(diff/gold_val > 0.5 || diff/acc_val > 0.5 || diff/gold_val < -0.5 || diff/acc_val < -0.5){
					if(j < x_dim)
						printf("ERROR: X Accelerator value: %f Golden value: %f index: %d iter: %d diff: %f\n", acc_val, gold_val, i * out_words_adj + j, i, diff);
					else
						printf("ERROR: P Accelerator value: %f Golden value: %f index: %d iter: %d diff: %f\n", acc_val, gold_val, i * out_words_adj + j, i, diff);
					errors++;
				}
			}
		}

	MSE /= ((x_dim + x_dim * x_dim) * iter);
	printf("Output MSE: %f \n", MSE);
	printf("Sum of errors is %u \n", errors);

	return errors;
}


/* User-defined code */
static void init_buffer(token_t *in, token_t * gold)
{
	int i;
	int j;

	/* for (i = 0; i < iter; i++) */
	/* 	for (j = 0; j < z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim; j++) */
	/* 		in[i * in_words_adj + j] = (token_t) j; */


	for(i = 0; i < iter; i++)
	{//z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim

		//Z
		for(j = 0; j < z_dim; j++)
		{
			if(i == 0)
				in[i * in_words_adj + j] = (token_t) measurements[NEURONS * (i+1) + j];
			else
				in[in_words_adj + (i-1) * in_words_adj_z + j] = (token_t) measurements[NEURONS * (i+1) + j];
			// if(i == 3)
			//     printf("Value of Z = %f index %d\n", measurements[NEURONS * (i+1) + j], in_words_adj + (i-1) * in_words_adj_z + j);
		}

		if(i == 0) //only for first iteration
		{
			//X
			for(; j < z_dim + x_dim; j++)
			{
				in[i * in_words_adj + j] = (token_t) initial[j - z_dim];
			}

			//P
			for(; j < z_dim + x_dim + x_dim * x_dim; j++)
			{
				in[i * in_words_adj + j] = (token_t) 0.0;
			}

			//F
			for(; j < z_dim + x_dim + x_dim * x_dim * 2; j++)
			{
				in[i * in_words_adj + j] = (token_t) A[j - (z_dim + x_dim + x_dim * x_dim)];
			}

			//Q
			for(; j < z_dim + x_dim + x_dim * x_dim * 3; j++)
			{
				in[i * in_words_adj + j] = (token_t) W[j - (z_dim + x_dim + x_dim * x_dim * 2)];
			}

			//R
			for(; j < z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim; j++)
			{
				in[i * in_words_adj + j] = (token_t) Q[j - (z_dim + x_dim + x_dim * x_dim * 3)];
			}

			//H
			for(; j < z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim; j++)
			{
				in[i * in_words_adj + j] = (token_t) H[j - (z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim)];
				//printf("Value of H = %f \n", measurements[NEURONS * (i+1) + j]);
			}
		}
	}


	/* for (i = 0; i < iter; i++) */
	/* 	for (j = 0; j < x_dim + x_dim * x_dim; j++) */
	/* 		gold[i * out_words_adj + j] = (token_t) j; */
	for(i = 0; i < iter; i++)
		for(j = 0; j < x_dim + x_dim * x_dim; j++)
		{
			if(j < x_dim)
				gold[i * out_words_adj + j] = (token_t) prediction[STATES * (i+1) + j];
			else
			{
				gold[i * out_words_adj + j] = (token_t) P_flat[i * x_dim * x_dim + (j - x_dim)];
			}
		}
}


/* User-defined code */
static void init_parameters()
{
	if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0) {
		in_words_adj = z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim;
		in_words_adj_z = z_dim;
		out_words_adj = x_dim + x_dim * x_dim;
	} else {
		in_words_adj = round_up(z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim, DMA_WORD_PER_BEAT(sizeof(token_t)));
		in_words_adj_z = round_up(z_dim, DMA_WORD_PER_BEAT(sizeof(token_t)));
		out_words_adj = round_up(x_dim + x_dim * x_dim, DMA_WORD_PER_BEAT(sizeof(token_t)));
	}
	in_len = in_words_adj + in_words_adj_z * (iter-1);
	out_len =  out_words_adj * (iter);
	in_size = in_len * sizeof(token_t);
	out_size = out_len * sizeof(token_t);
	out_offset = in_len;
	size = (out_offset * sizeof(token_t)) + out_size;
}


int main(int argc, char **argv)
{
	int errors;

	token_t *gold;
	token_t *buf;

	init_parameters();

	buf = (token_t *) esp_alloc(size);
	cfg_000[0].hw_buf = buf;

	gold = malloc(out_size);

	init_buffer(buf, gold);

	printf("\n====== %s ======\n\n", cfg_000[0].devname);
	/* <<--print-params-->> */
	printf("  .iter = %d\n", iter);
	printf("  .x_dim = %d\n", x_dim);
	printf("  .z_dim = %d\n", z_dim);
	printf("\n  ** START **\n");

	esp_run(cfg_000, NACC);

	printf("\n  ** DONE **\n");

	errors = validate_buffer(&buf[out_offset], gold);

	free(gold);
	esp_free(buf);

	if (!errors)
		printf("+ Test PASSED\n");
	else
		printf("+ Test FAILED\n");

	printf("\n====== %s ======\n\n", cfg_000[0].devname);

	return errors;
}

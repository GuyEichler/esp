// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "libesp.h"
#include "cfg_flex_soma.h"

#include "A_array_soma.h"
#include "H_array_soma.h"
#include "W_array_soma.h"
#include "Q_array_soma.h"
#include "initial_state_array_soma.h"
#include "measurements_array_soma.h"
#include "prediction_array_soma.h"
#include "real_array_soma.h"
#include "P_array_soma.h"

#include "SSK_gain.h"

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
	float MAE = 0.0;
	float max_diff = 0.0;

	for (i = 0; i < iter; i++)
		for (j = 0; j < (x_dim + x_dim * x_dim)*chunks; j++)
		{

			token_t gold_val = gold[i * out_words_adj + j];
			token_t acc_val = out[i * out_words_adj + j];

			token_t diff;
			if(gold_val < acc_val)
				diff = acc_val - gold_val;
			else
				diff = gold_val - acc_val;

			MSE += diff * diff;
			MAE += diff;

			token_t norm_diff = diff/acc_val > diff/gold_val ? diff/acc_val : diff/gold_val;
			if(norm_diff > max_diff)
				max_diff = norm_diff;

			/* if(j%(x_dim + x_dim*x_dim) < x_dim) */
			/* 	printf("NO ERROR: X Accelerator value: %f Golden value: %f index: %d iter: %d diff: %f\n", acc_val, gold_val, i * out_words_adj + j, i, diff); */
			/* else */
			/* 	printf("N ERROR: P Accelerator value: %f Golden value: %f index: %d iter: %d diff: %f\n", acc_val, gold_val, i * out_words_adj + j, i, diff); */

			if (gold[i * out_words_adj + j] != out[i * out_words_adj + j])
			{
				if(diff/gold_val > 0.1 || diff/acc_val > 0.1 || diff/gold_val < -0.1 || diff/acc_val < -0.1){
					/* if(j < x_dim) */
					/* 	printf("ERROR: X Accelerator value: %f Golden value: %f index: %d iter: %d diff: %f\n", acc_val, gold_val, i * out_words_adj + j, i, diff); */
					/* else */
					/* 	printf("ERROR: P Accelerator value: %f Golden value: %f index: %d iter: %d diff: %f\n", acc_val, gold_val, i * out_words_adj + j, i, diff); */
					errors++;
				}
			}
		}

	MSE /= ((x_dim + x_dim * x_dim) * iter*chunks);
	MAE /= ((x_dim + x_dim * x_dim) * iter*chunks);
	printf("Output MSE: %.20f \n", MSE);
	printf("Output MAE: %.20f \n", MAE);
	printf("Maximum difference: %f percent \n", max_diff);
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

		j = 0;

		if(i == 0) //only for first iteration
		{
			//X
			for(; j < x_dim; j++)
			{
				in[i * in_words_adj + j] = (token_t) initial[j];
			}

			//P
			for(; j < x_dim + x_dim * x_dim; j++)
			{
				in[i * in_words_adj + j] = (token_t) 0.0;
			}

			//F
			for(; j < x_dim + x_dim * x_dim * 2; j++)
			{
				in[i * in_words_adj + j] = (token_t) A[j - (x_dim + x_dim * x_dim)];
			}

			//Q
			for(; j < x_dim + x_dim * x_dim * 3; j++)
			{
				in[i * in_words_adj + j] = (token_t) W[j - (x_dim + x_dim * x_dim * 2)];
			}

			/* //R */
			/* for(; j < x_dim + x_dim * x_dim * 3 + z_dim * z_dim; j++) */
			/* { */
			/* 	in[i * in_words_adj + j] = (token_t) Q[j - (x_dim + x_dim * x_dim * 3)]; */
			/* } */

			//H
			for(; j < x_dim + x_dim * x_dim * 3 + z_dim * x_dim; j++)
			{
				in[i * in_words_adj + j] = (token_t) H[j - (x_dim + x_dim * x_dim * 3)];
				//printf("Value of H = %f \n", measurements[NEURONS * (i+1) + j]);
			}

			//SSK
			for(; j < x_dim + x_dim * x_dim * 3 + 2 * z_dim * x_dim; j++)
			{
				in[i * in_words_adj + j] = (token_t) SSK_gain[j - (x_dim + x_dim * x_dim * 3 + z_dim * x_dim)];// + rand_val;
				//printf("Value of H = %f \n", measurements[NEURONS * (i+1) + j]);
			}
		}

		unsigned base_index = (x_dim + x_dim * x_dim * 3 + 2 * z_dim * x_dim);

		//Z
		if(i == 0)
			for(; j < base_index + z_dim * chunks; j++)
			{

				in[j] = (token_t) measurements[z_dim * (i+1) + j - base_index];

				/* in[in_words_adj + (i-1) * in_words_adj_z + j] = (token_t) measurements[NEURONS * (i+1) + j]; */
			}
		else
			for(; j < z_dim * chunks; j++)
			{
				in[in_words_adj + i * in_words_adj_z + j] = (token_t) measurements[z_dim * (i * chunks + 1) + j];
				//int32_t val = float_to_fixed32(measurements[NEURONS * i + j], 3);
			// if(i == 3)
			//printf("Value of Z = %d index %d \n", val, i * in_words_adj + j);
			}
	}


	/* for (i = 0; i < iter; i++) */
	/* 	for (j = 0; j < x_dim + x_dim * x_dim; j++) */
	/* 		gold[i * out_words_adj + j] = (token_t) j; */
	for(i = 0; i < iter*chunks; i++)
		for(j = 0; j < x_dim + x_dim * x_dim; j++)
		{
			if(j < x_dim)
				gold[i * out_words_adj/chunks + j] = (token_t) prediction[x_dim * (i+1) + j];
			else
			{
				gold[i * out_words_adj/chunks + j] = (token_t) P_flat[i * x_dim * x_dim + (j - x_dim)];
			}
		}
}


/* User-defined code */
static void init_parameters()
{
	if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0) {
		in_words_adj = x_dim + x_dim * x_dim * 3 + 2 * z_dim * x_dim;
		in_words_adj_z = z_dim * chunks;
		out_words_adj = (x_dim + x_dim * x_dim) * chunks;
	} else {
		in_words_adj = round_up(x_dim + x_dim * x_dim * 3 + 2 * z_dim * x_dim, DMA_WORD_PER_BEAT(sizeof(token_t)));
		in_words_adj_z = round_up(z_dim * chunks, DMA_WORD_PER_BEAT(sizeof(token_t)));
		out_words_adj = round_up((x_dim + x_dim * x_dim) * chunks, DMA_WORD_PER_BEAT(sizeof(token_t)));
	}
	in_len = in_words_adj + in_words_adj_z * iter;
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

	if(argc < 7) {
		printf("running x_dim=%d z_dim=%d iter=%d chunks=%d inv_num=%d inv_reset=%d \n", x_dim, z_dim, iter, chunks, inv_num, inv_reset);
	}
	else {
		const char* x_str = argv[1];
		x_dim = atoi(x_str);
		((struct kalman_sskg_vivado_access*) cfg_000[0].esp_desc)->x_dim = x_dim;
		const char* z_str = argv[2];
		z_dim = atoi(z_str);
		((struct kalman_sskg_vivado_access*) cfg_000[0].esp_desc)->z_dim = z_dim;
		const char* iter_str = argv[3];
		iter = atoi(iter_str);
		((struct kalman_sskg_vivado_access*) cfg_000[0].esp_desc)->iter = iter;
		const char* chunks_str = argv[4];
		chunks = atoi(chunks_str);
		((struct kalman_sskg_vivado_access*) cfg_000[0].esp_desc)->chunks = chunks;
		const char* inv_num_str = argv[5];
		inv_num = atoi(inv_num_str);
		((struct kalman_sskg_vivado_access*) cfg_000[0].esp_desc)->inv_num = inv_num;
		const char* inv_reset_str = argv[6];
		inv_reset = atoi(inv_reset_str);
		((struct kalman_sskg_vivado_access*) cfg_000[0].esp_desc)->inv_reset = inv_reset;
		printf("running x_dim=%d z_dim=%d iter=%d chunks=%d inv_num=%d inv_reset=%d \n", x_dim, z_dim, iter, chunks, inv_num, inv_reset);
	}


	init_parameters();

	buf = (token_t *) esp_alloc(size);
	cfg_000[0].hw_buf = buf;

	gold = malloc(out_size);

	init_buffer(buf, gold);

	printf("\n====== %s ======\n\n", cfg_000[0].devname);
	/* <<--print-params-->> */
	printf("  .inv_reset = %d\n", ((struct kalman_sskg_vivado_access*) cfg_000[0].esp_desc)->inv_reset);
	printf("  .inv_num = %d\n", ((struct kalman_sskg_vivado_access*) cfg_000[0].esp_desc)->inv_num);
	printf("  .chunks = %d\n", ((struct kalman_sskg_vivado_access*) cfg_000[0].esp_desc)->chunks);
	printf("  .iter = %d\n", ((struct kalman_sskg_vivado_access*) cfg_000[0].esp_desc)->iter);
	printf("  .x_dim = %d\n", ((struct kalman_sskg_vivado_access*) cfg_000[0].esp_desc)->x_dim);
	printf("  .z_dim = %d\n", ((struct kalman_sskg_vivado_access*) cfg_000[0].esp_desc)->z_dim);
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

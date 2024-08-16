// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "libesp.h"
#include "cfg_small.h"

/* #include "A_array.h" */
/* #include "H_array.h" */
/* #include "W_array.h" */
/* #include "Q_array.h" */
/* #include "initial_state_array.h" */
/* #include "measurements_array.h" */
/* #include "prediction_array.h" */
/* #include "real_array.h" */
/* #include "P_array.h" */
#include "small.h"

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

			if(j%(x_dim + x_dim*x_dim) < x_dim)
				printf("NO ERROR: X Accelerator value: %f Golden value: %f index: %d iter: %d diff: %f\n", acc_val, gold_val, i * out_words_adj + j, i, diff);
			else
				printf("NO ERROR: P Accelerator value: %f Golden value: %f index: %d iter: %d diff: %f\n", acc_val, gold_val, i * out_words_adj + j, i, diff);

			if (gold[i * out_words_adj + j] != out[i * out_words_adj + j])
			{
				if(diff / gold_val > 0.1 || diff / acc_val > 0.1 || diff / gold_val < -0.1 || diff / acc_val < -0.1){
					if(j < x_dim)
						printf("ERROR: X Accelerator value: %f Golden value: %f index: %d iter: %d diff: %f\n", acc_val, gold_val, i * out_words_adj + j, i, diff);
					else
						printf("ERROR: P Accelerator value: %f Golden value: %f index: %d iter: %d diff: %f\n", acc_val, gold_val, i * out_words_adj + j, i, diff);
					errors++;
				}
			}
		}

	MSE /= ((x_dim + x_dim * x_dim) * TIME_STAMPS);
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

	/* printf("Initializing buffer...\n"); */

	for(i = 0; i < iter; i++)
	{//z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim

		/* printf("Initializing model...\n"); */

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
		/* printf("Before loading Z \n"); */
		//Z
		if(i == 0)
			for(; j < base_index + z_dim * chunks; j++)
			{
				in[j] = (token_t) measurements[NEURONS * i + j - base_index];

				/* printf("Value of Z = %f index %d \n", in[j], j); */
			}
		else
			for(; j < z_dim * chunks; j++)
			{
				in[in_words_adj + i * in_words_adj_z + j] = (token_t) measurements[NEURONS * i * chunks + j];

				/* printf("Value of Z = %f index %d \n", in[in_words_adj + i * in_words_adj_z + j], in_words_adj + i * in_words_adj_z + j); */
			}
		/* printf("After loading Z \n"); */
	}


	/* for (i = 0; i < iter; i++) */
	/* 	for (j = 0; j < x_dim + x_dim * x_dim; j++) */
	/* 		gold[i * out_words_adj + j] = (token_t) j; */
	for(i = 0; i < TIME_STAMPS; i++)
		for(j = 0; j < x_dim + x_dim * x_dim; j++)
		{
			if(j < x_dim)
				gold[i * out_words_adj/chunks + j] = (token_t) prediction[STATES * i + j];
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

	init_parameters();

	buf = (token_t *) esp_alloc(size);
	if(buf == NULL) printf("Mem error!\n");

	cfg_000[0].hw_buf = buf;

	gold = malloc(out_size);

	init_buffer(buf, gold);

	printf("\n====== %s ======\n\n", cfg_000[0].devname);
	/* <<--print-params-->> */
	printf("  .inv_reset = %d\n", ((struct kalman_sskg_vivado_access *)cfg_000[0].esp_desc)->inv_reset);
	printf("  .inv_num = %d\n", ((struct kalman_sskg_vivado_access *)cfg_000[0].esp_desc)->inv_num);
	printf("  .chunks = %d\n", ((struct kalman_sskg_vivado_access *)cfg_000[0].esp_desc)->chunks);
	printf("  .iter = %d\n", ((struct kalman_sskg_vivado_access *)cfg_000[0].esp_desc)->iter);
	printf("  .x_dim = %d\n", ((struct kalman_sskg_vivado_access *)cfg_000[0].esp_desc)->x_dim);
	printf("  .z_dim = %d\n", ((struct kalman_sskg_vivado_access *)cfg_000[0].esp_desc)->z_dim);
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

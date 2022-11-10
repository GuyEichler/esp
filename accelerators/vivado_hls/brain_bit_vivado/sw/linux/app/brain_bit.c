// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "libesp.h"
#include "cfg.h"
#include "input_full.h"

static unsigned in_words_adj;
static unsigned out_words_adj;
static unsigned in_len;
static unsigned out_len;
static unsigned in_size;
static unsigned out_size;
static unsigned out_offset;
static unsigned size;
static int key_counter;

/* User-defined code */
static int validate_buffer(token_t *out, token_t *gold)
{
	int i;
	int j;
	unsigned errors = 0;

	int skip = 0;
	key_counter = 0;

	for (i = 0; i < key_batch; i++)
		for (j = 0; j < key_length; j++){
			if(key_counter == key_num) break;
			unsigned index = i * out_words_adj + j;
			/* token_t val = out[index - skip]; */
                        int bit = (index - skip) % 32;
                        int word = (index - skip) >> 5;
                        token_t mask = out[word] >> bit;
                        /* printf("out val is %x for word %d bit %d\n", out[word], word, bit); */
                        token_t val = mask & 1;
			token_t gold_val = gold[index];
			unsigned reduce = (ceil((float)skip/key_length));
			if(gold_val != 3){
				if(!(i == key_batch - reduce && (j > skip - 1) )){
					printf("Calculated value %x Golden value %d for index %d \n",
						val, gold_val, (index-skip));
					if (gold_val != val){
						errors++;
						printf("ERROR\n");
					}
				}
			}
			else{
				printf("SKIPPING\n");
				skip += 1;
			}

			if((index - skip + 1) % key_length == 0 && index != 0){
				key_counter++;
				printf("\n----------KEY %d DONE----------\n", key_counter);
                                printf("\nKEY IS: [ ");
                                for(int k = key_length / 32 - 1; k >= 0; k--)
					printf("0x%x ", out[word-k]);
                                printf("]\n\n");
			}
		}

	return errors;
}


/* User-defined code */
static void init_buffer(token_t *in, token_t * gold)
{
	int i;
	int j;

	for (i = 0; i < key_batch; i++)
		for (j = 0; j < key_length; j++){
                        float val = val_arr[i * in_words_adj + j];
                        in[i * in_words_adj + j] = (token_t) float_to_fixed32(val, 12);
                        //in[i * in_words_adj + j] = (token_t) val;
                        /* printf("Generated value %f\n", fixed32_to_float(in[i * out_words_adj + j] , 12)); */
                }

        for (i = 0; i < key_batch; i++)
		for (j = 0; j < key_length; j++){
			float val = val_arr[i * in_words_adj + j];
                        bool filter = (fabs((float)val - avg) >= Rs);
			if(!filter){
				int32_t result = floor((float)(((val - (avg - Rs)) / (2*Rs)) * L));
				result = result % 2;
				gold[i * out_words_adj + j] = (token_t) result;
                                //printf("Generated golden value %d\n", gold[i * out_words_adj + j]);
			}
			else{
				gold[i * out_words_adj + j] = 3;
			}
		}
}


/* User-defined code */
static void init_parameters()
{
	if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0) {
		in_words_adj = key_length;
		out_words_adj = key_length;
	} else {
		in_words_adj = round_up(key_length, DMA_WORD_PER_BEAT(sizeof(token_t)));
		out_words_adj = round_up(key_length, DMA_WORD_PER_BEAT(sizeof(token_t)));
	}
	in_len = in_words_adj * (key_batch);
	out_len =  out_words_adj * (key_batch);
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
	printf("  .avg = %f\n", avg);
	printf("  .key_length = %d\n", key_length);
	printf("  .std = %f\n", std);
	printf("  .R = %f\n", R);
	printf("  .L = %d\n", L);
	printf("  .key_batch = %d\n", key_batch);
	printf("  .key_num = %d\n", key_num);
	printf("\n  ** START **\n");

        unsigned avg_u = *avg_ptr;
        unsigned std_u = *std_ptr;
        unsigned R_u = *R_ptr;

        ((struct brain_bit_vivado_access*) cfg_000[0].esp_desc)->avg = avg_u;
        ((struct brain_bit_vivado_access*) cfg_000[0].esp_desc)->std = std_u;
        ((struct brain_bit_vivado_access*) cfg_000[0].esp_desc)->R = R_u;

        /* printf("Value avg: %x\n", ((struct brain_bit_vivado_access*) cfg_000[0].esp_desc)->avg); */
        /* printf("Value std: %x\n", ((struct brain_bit_vivado_access*) cfg_000[0].esp_desc)->std); */
        /* printf("Value R: %x\n", ((struct brain_bit_vivado_access*) cfg_000[0].esp_desc)->R); */

	/* token_t* out_location = &buf[out_offset]; */
	/* for(int i = 0; i < out_len; i++) */
	/* 	out_location[i] = 0; */

	esp_run(cfg_000, NACC);

	printf("\n  ** DONE **\n");

	errors = validate_buffer(&buf[out_offset], gold);

	free(gold);
	esp_free(buf);

        float total = 100 * (float) errors / (key_length*key_batch);
	if (total <= 1){
                printf("- Keys generated: %d\n", key_counter);
		printf("- Number of bit errors: %d\n", errors);
		printf("+ Test PASSED\n");
	}
	else{
                printf("+ Test FAILED with %f%% errors\n", total);
	}

	printf("\n====== %s ======\n\n", cfg_000[0].devname);

	return errors;
}

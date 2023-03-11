// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0

#include "libesp.h"
#include "cfg.h"
#include <stdio.h>
#include <sys/time.h>
#include "cfg_rand.h"
#include "input_1mil_full.h"

typedef int32_t token_t;

//#define SYS_RAND
#define PRINT 0
#define DEV_RAND
#define KEY_LEN 256000
#define KEY_BYTES 32
#define MAX_ITER  1

#define RAND_MODE SYS_RAND

const unsigned num_rand_words = KEY_LEN/ (sizeof(token_t) * 8);
/* const unsigned num_rand_words = KEY_LEN/ (sizeof(token_t)); */
const unsigned key_len_bytes = KEY_LEN / 8;
unsigned char *keys;
char *key_chars;
struct timeval hw_begin_0, hw_end_0;
int i, k;
long seconds = 0, microseconds = 0;
double elapsed = 0.0, avg_elapsed_us = 0.0, total_time = 0.0;

//brain_bit
static unsigned in_words_adj;
static unsigned out_words_adj;
static unsigned in_len;
static unsigned out_len;
static unsigned in_size;
static unsigned out_size;
static unsigned out_offset;
static unsigned size;
static int key_counter;

static int validate_buffer(token_t *out, token_t *gold)
{
	int i;
	int j;
	unsigned errors = 0;

	int skip = 0;
	key_counter = 0;
	/* int val_counter = 0; */
	int offset = 0;
	bool done = false;

	for (i = 0; i < key_batch; i++)
		for (j = 0; j < key_length; j++){
                        if(key_counter != key_num){
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
					if (gold_val != val){
                                                printf("Calculated value %x Golden value %d for index %d \n",val, gold_val, (index-skip));
						errors++;
						printf("ERROR\n");
					}
				}
			}
			else{
				/* printf("SKIPPING\n"); */
				skip += 1;
			}

			if((index - skip + 1) % (key_length*(key_counter+1)) == 0 && index != 0){
				key_counter++;
				/* printf("\n----------KEY %d DONE----------\n", key_counter); */
				/* printf("\nKEY IS: [ "); */
				/* for(int k = key_length / 32 - 1; k >= 0; k--) */
				/* 	printf("0x%x ", out[word-k]); */
				/* printf("]\n\n"); */
			}
                        }
			else if(!done){
				done = true;
				unsigned max_chunk = 1024;
				offset = i * out_words_adj + j - (skip % key_length);
				if(key_length > max_chunk && skip > key_length)
					offset += key_length - max_chunk;
				else if(key_length > max_chunk && skip < key_length)
					offset -= key_length % max_chunk;
			}
		}

	int index_offset = key_num * key_length / DATA_BITWIDTH;

	for(int i = 0; i < val_num; i++){
		unsigned index = index_offset + i + DATA_BITWIDTH / key_length;
                token_t val = out[index];
                token_t gold_val = (token_t) float_to_fixed32(val_arr[offset+i], 12);

                if(val != gold_val){
			printf("Calculated value %x Golden value %x for index %d \n", val, gold_val, index);
			printf("ERROR\n");
			errors += key_length;
		}
                /* if((index + 1) % key_length == 0 && index != 0) */
                /*         val_counter++; */
         }


	return errors;
}

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

int run_dev_rand()
{
    FILE *f;

    keys = malloc(key_len_bytes);
    if(keys == NULL) {
        printf("Error allocating mem \n");
        return -1;
    }

    //average time to read random numbers
    for(k = 0; k < MAX_ITER; k++) {

        //start timer
        gettimeofday(&hw_begin_0, 0);

        f = fopen("/dev/random", "r");
        if(f == NULL) {
            printf("Error opening /dev/random \n");
            return -1;
        }

        fread((void *) keys, key_len_bytes, 1, f);

        fclose(f);

        //end timer
        gettimeofday(&hw_end_0, 0);

        //print rand number is hex
        for(i = 0; i < num_rand_words; i++)
	    if(PRINT)
	        printf("%x", (unsigned int) keys[i * sizeof(token_t)]);

        printf("\n");

        seconds = hw_end_0.tv_sec - hw_begin_0.tv_sec;
        microseconds = hw_end_0.tv_usec - hw_begin_0.tv_usec;
        elapsed += seconds + microseconds * 1e-6;
    }

    free(keys);

    total_time = (double) elapsed * 1e6;
    avg_elapsed_us = total_time / MAX_ITER;

    printf("Randomness source: /dev/random\n");
    printf("Number of random bits: %d\n", KEY_LEN*MAX_ITER);
    printf("Time measured: %u Iterations, avg time %.4f us.\n", MAX_ITER, avg_elapsed_us);

    return 0;
}

int run_dev_urand()
{
    FILE *f;

    keys = malloc(key_len_bytes);
    if(keys == NULL) {
        printf("Error allocating mem \n");
        return -1;
    }

    //average time to read random numbers
    for(k = 0; k < MAX_ITER; k++) {

        //start timer
        gettimeofday(&hw_begin_0, 0);

        f = fopen("/dev/urandom", "r");
        if(f == NULL) {
            printf("Error opening /dev/urandom \n");
            return -1;
        }

        fread((void *) keys, key_len_bytes, 1, f);

        fclose(f);

        //end timer
        gettimeofday(&hw_end_0, 0);

        //print rand number is hex
        for(i = 0; i < num_rand_words; i++)
	    if(PRINT)
	        printf("%x", (unsigned int) keys[i * sizeof(token_t)]);

        printf("\n");

        seconds = hw_end_0.tv_sec - hw_begin_0.tv_sec;
        microseconds = hw_end_0.tv_usec - hw_begin_0.tv_usec;
        elapsed += seconds + microseconds * 1e-6;
    }

    free(keys);

    total_time = (double) elapsed * 1e6;
    avg_elapsed_us = total_time / MAX_ITER;

    printf("Randomness source: /dev/urandom\n");
    printf("Number of random bits: %d\n", KEY_LEN*MAX_ITER);
    printf("Time measured: %u Iterations, avg time %.4f us.\n", MAX_ITER, avg_elapsed_us);

    return 0;
}

int run_rand()
{
   keys = malloc(key_len_bytes);

    for (k = 0; k < MAX_ITER; k++) {
        gettimeofday(&hw_begin_0, 0);

        srand(k);
        for(i = 0; i < num_rand_words; i++){
            keys[i * sizeof(token_t)] = rand();
        }

        gettimeofday(&hw_end_0, 0);

        for(i = 0; i < num_rand_words; i++)
	    if(PRINT)
	        printf("%x", (unsigned int) keys[i * sizeof(token_t)]);

        printf("\n");
        seconds = hw_end_0.tv_sec - hw_begin_0.tv_sec;
        microseconds = hw_end_0.tv_usec - hw_begin_0.tv_usec;
        elapsed += seconds + microseconds*1e-6;
    }

    avg_elapsed_us = (double) (elapsed * 1e6) / MAX_ITER;

    printf("Randomness source: rand()\n");
    printf("Number of random bits: %d\n", KEY_LEN*MAX_ITER);
    printf("Time measured: %u Iterations, avg time %.4f us.\n", MAX_ITER, avg_elapsed_us);

    free(keys);
    return 0;
}

int run_brain()
{

	token_t *gold;
	token_t *buf;

	init_parameters();
	buf = (token_t *) esp_alloc(size);
	cfg_000[0].hw_buf = buf;
	gold = malloc(out_size);

	init_buffer(buf, gold);

	//Set up key size (number of random bits) and number of keys to generate (iterations)
	/* key_length = KEY_LEN; */
	/* key_num = MAX_ITER; */

	printf("\n====== %s ======\n\n", cfg_000[0].devname);
	/* <<--print-params-->> */
	printf("  .avg = %f\n", avg);
	printf("  .key_length = %d\n", KEY_LEN);
	printf("  .std = %f\n", std);
	printf("  .R = %f\n", R);
	printf("  .L = %d\n", L);
	printf("  .key_batch = %d\n", key_batch);
	printf("  .key_num = %d\n", MAX_ITER);
	printf("\n  ** START **\n");

        unsigned avg_u = *avg_ptr;
        unsigned std_u = *std_ptr;
        unsigned R_u = *R_ptr;

        ((struct brain_bit_vivado_access*) cfg_000[0].esp_desc)->avg = avg_u;
        ((struct brain_bit_vivado_access*) cfg_000[0].esp_desc)->std = std_u;
        ((struct brain_bit_vivado_access*) cfg_000[0].esp_desc)->R = R_u;
        ((struct brain_bit_vivado_access*) cfg_000[0].esp_desc)->key_num = MAX_ITER;
        ((struct brain_bit_vivado_access*) cfg_000[0].esp_desc)->key_length = KEY_LEN;

	token_t* out_location = &buf[out_offset];
	for(int i = 0; i < out_len; i++)
		out_location[i] = 0;

	unsigned long long total_time = 0;
	/* for(int k = 0; k < N_runs; k++){ */
	esp_run_no_print(cfg_000, NACC);
	total_time = cfg_000[0].hw_ns;
	/* } */

	int errors = 0;
	if(PRINT)
            errors = validate_buffer(&buf[out_offset], gold);


	avg_elapsed_us = ((double)total_time / 1e3) / MAX_ITER;

	printf("Randomness source: brain_bit\n");
	printf("Number of random bits: %d\n", KEY_LEN*MAX_ITER);
	printf("Time measured: %u Iterations, avg time %.4f us.\n", MAX_ITER, avg_elapsed_us);

	free(keys);
	return errors;
}

int main(int argc, char *argv[])
{
    if(argc < 2 || *argv[1] == '1') {
        printf("running /dev/urandom \n");
        run_dev_urand();
    }
    else if (*argv[1] == '2') {
        printf("running /dev/random \n");
        run_dev_rand();
    }
    else if (*argv[1] == '3') {
        printf("running rand \n");
        run_rand();
    }
    else if (*argv[1] == '4') {
        printf("running brain_bit \n");
        run_brain();
    }

    return 0;
}

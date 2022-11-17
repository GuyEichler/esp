// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "libesp.h"
#include "cfg.h"
#include "brain_bit_input_full.h"
#include "aes_data.h"

#include <sys/time.h>

// brain_bit:
static unsigned in_words_adj;
static unsigned out_words_adj;
static unsigned in_len;
static unsigned out_len;
static unsigned in_size;
static unsigned out_size;
static unsigned out_offset;
static unsigned size;
static int key_counter;

// aes:
#define N_BATCH 9
// static unsigned in_words_adj;
// static unsigned out_words_adj;
// static unsigned in_len;
// static unsigned out_len;
// static unsigned in_size;
// static unsigned out_size;
// static unsigned out_offset;
static unsigned size_bytes;

// static unsigned tag_bytes;
/* static unsigned aad_bytes; */
/* static unsigned in_bytes; */
/* static unsigned out_bytes; */
/* static unsigned iv_bytes; */
/* static unsigned key_bytes; */
/* static unsigned encryption; */
/* static unsigned oper_mode; */

/* static unsigned key_words; */
/* static unsigned iv_words; */
static unsigned in_words;
static unsigned out_words;
/* static unsigned aad_words; */
/* static unsigned tag_words; */

/* static unsigned key_size; */
/* static unsigned iv_size; */
// static unsigned in_size;
// static unsigned out_size;
/* static unsigned aad_size; */
/* static unsigned tag_size; */

// input_size for each iteration
const unsigned array_size[12] = {16, 32, 48, 64, 80, 96, 112, 128, 144};

static void init_parameters_brain_bit()
{
	if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0)
	{
		in_words_adj = key_length;
		out_words_adj = key_length;
	}
	else
	{
		in_words_adj = round_up(key_length, DMA_WORD_PER_BEAT(sizeof(token_t)));
		out_words_adj = round_up(key_length, DMA_WORD_PER_BEAT(sizeof(token_t)));
	}
	in_len = in_words_adj * (key_batch);
	out_len = out_words_adj * (key_batch);
	in_size = in_len * sizeof(token_t);
	out_size = out_len * sizeof(token_t);
	out_offset = in_len;
	size = (out_offset * sizeof(token_t)) + out_size;
}

static void init_buffer_brain_bit(token_t *in, token_t *gold)
{
	int i;
	int j;

	for (i = 0; i < key_batch; i++)
		for (j = 0; j < key_length; j++)
		{
			float val = val_arr[i * in_words_adj + j];
			in[i * in_words_adj + j] = (token_t)float_to_fixed32(val, 12);
			// in[i * in_words_adj + j] = (token_t) val;
			/* printf("Generated value %f\n", fixed32_to_float(in[i * out_words_adj + j] , 12)); */
		}

	for (i = 0; i < key_batch; i++)
		for (j = 0; j < key_length; j++)
		{
			float val = val_arr[i * in_words_adj + j];
			bool filter = (fabs((float)val - avg) >= Rs);
			if (!filter)
			{
				int32_t result = floor((float)(((val - (avg - Rs)) / (2 * Rs)) * L));
				result = result % 2;
				gold[i * out_words_adj + j] = (token_t)result;
				// printf("Generated golden value %d\n", gold[i * out_words_adj + j]);
			}
			else
			{
				gold[i * out_words_adj + j] = 3;
			}
		}
}

static int validate_buffer_brain_bit(token_t *out, token_t *gold)
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
	{
		for (j = 0; j < key_length; j++)
		{
			if (key_counter != key_num)
			{
				unsigned index = i * out_words_adj + j;
				/* token_t val = out[index - skip]; */
				int bit = (index - skip) % 32;
				int word = (index - skip) >> 5;
				token_t mask = out[word] >> bit;
				/* printf("out val is %x for word %d bit %d\n", out[word], word, bit); */
				token_t val = mask & 1;
				token_t gold_val = gold[index];
				unsigned reduce = (ceil((float)skip / key_length));
				if (gold_val != 3)
				{
					if (!(i == key_batch - reduce && (j > skip - 1)))
					{
						if (gold_val != val)
						{
							printf("Calculated value %x Golden value %d for index %d \n", val, gold_val, (index - skip));
							errors++;
							printf("ERROR\n");
						}
					}
				}
				else
				{
					printf("SKIPPING\n");
					skip += 1;
				}

				if ((index - skip + 1) % key_length == 0 && index != 0)
				{
					key_counter++;
					printf("----------KEY %d DONE----------\n", key_counter);
					printf("KEY IS: [ ");
					for (int k = key_length / 32 - 1; k >= 0; k--)
						printf("0x%x ", out[word - k]);
					printf("]\n");
				}
			}
			else if (!done)
			{
				done = true;
				offset = i * out_words_adj + j - skip;
			}
		}
	}

	int index_offset = key_num * key_length / DATA_BITWIDTH;

	for (int i = 0; i < val_num; i++)
	{
		unsigned index = index_offset + i + DATA_BITWIDTH / key_length;
		token_t val = out[index];
		token_t gold_val = (token_t)float_to_fixed32(val_arr[offset + i], 12);

		if (val != gold_val)
		{
			printf("Calculated value %x Golden value %x for index %d \n", val, gold_val, index);
			printf("ERROR\n");
			errors += key_length;
		}
		/* if((index + 1) % key_length == 0 && index != 0) */
		/*         val_counter++; */
	}

	return errors;
}

static void init_parameters_aes(unsigned indx)
{
	in_words = ecb_raw_encrypt_plaintext_words[indx] +
			   ecb_raw_encrypt_key_words[indx];
	out_words = ecb_raw_encrypt_plaintext_words[indx];

	if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0)
	{
		in_words_adj = in_words;
		out_words_adj = out_words;
	}
	else
	{
		in_words_adj = round_up(in_words, DMA_WORD_PER_BEAT(sizeof(token_t)));
		out_words_adj = round_up(out_words, DMA_WORD_PER_BEAT(sizeof(token_t)));
	}

	in_len = in_words_adj;
	out_len = out_words_adj;
	in_size = in_len * sizeof(token_t);
	out_size = out_len * sizeof(token_t);
	out_offset = in_len;
	size_bytes = (out_offset * sizeof(token_t)) + out_size;

	printf("%s: in_words_adj  = %u\n", __func__, in_words_adj);
	printf("%s: out_words_adj = %u\n", __func__, out_words_adj);
	printf("%s: in_len        = %u\n", __func__, in_len);
	printf("%s: out_len       = %u\n", __func__, out_len);
	printf("%s: in_size       = %u\n", __func__, in_size);
	printf("%s: out_size      = %u\n", __func__, out_size);
	printf("%s: out_offset    = %u\n", __func__, out_offset);
	printf("%s: size_bytes    = %u\n", __func__, size_bytes);
}

static void init_buffer_aes(token_t *in, token_t *gold, token_t *out, unsigned indx)
{
	int i;
	int j;

	printf("init_buffer_aes: raw_encrypt_plaintext_words     [%u] %u\n", indx, ecb_raw_encrypt_plaintext_words[indx]);
	printf("init_buffer_aes: ecb_raw_encrypt_ciphertext_words[%u] %u\n", indx, ecb_raw_encrypt_ciphertext_words[indx]);

	for (j = 0; j < ecb_raw_encrypt_key_words[indx]; j++)
	{
		in[j] = ecb_raw_encrypt_key[indx][j];
		// printf("INFO: raw_encrypt_key[%u][%u] | %x\n", indx, j, in[j]);
	}

	for (i = 0; i < ecb_raw_encrypt_plaintext_words[indx]; i++, j++)
	{
		in[j] = ecb_raw_encrypt_plaintext[indx][i];
		// printf("INFO: raw_encrypt_plaintext[%u][%u] | inputs[%u]@%p %x\n", indx, i, j, in + j, in[j]);
	}

	for (j = 0; j < ecb_raw_encrypt_ciphertext_words[indx]; j++)
	{
		gold[j] = ecb_raw_encrypt_ciphertext[indx][j];
		// printf("INFO: raw_encrypt_ciphertext[%u][%u] %x\n", indx, j, gold[j]);
	}
}

static int validate_buffer_aes(token_t *out, token_t *gold, unsigned indx)
{
	/* int i; */
	int j;
	unsigned errors = 0;

	printf("INFO:  gold output data @%p\n", gold);
	printf("INFO:       output data @%p\n", out);

	for (j = 0; j < ecb_raw_encrypt_ciphertext_words[indx] /* + 0x10 */; j++)
	{
		token_t gold_data = gold[j];
		token_t out_data = out[j];

		if (out_data != gold_data)
		{
			errors++;
			printf("INFO: [%u] @%p %x (%x) %s\n", j, out + j, out_data, gold_data, ((out_data != gold_data) ? " !!!" : ""));
		}
	}

	printf("INFO: total errors %u\n", errors);

	return errors;
}

static void set_aes_in_from_brain_bit_out(token_t *in_aes, token_t *out_brain)
{
	int i;
	int j;

	int key_num_words = key_num / N_TESTS;
	int val_num_words = val_num / N_TESTS;

	for (i = 0; i < key_num_words; i++)
	{
		in_aes[i] = out_brain[i];
	}
	for (j = 0; j < val_num_words; j++)
	{
		in_aes[i + j] = out_brain[j + i];
	}
}

static void init_buffer_aes_from_brain(token_t *in, token_t *aes_key, token_t *aes_val, token_t *out, unsigned indx)
{

	int i;
	int j;
	int key_words = key_length / DATA_BITWIDTH;
	int val_words = val_num / N_TESTS;

	for (j = 0; j < key_words; j++)
	{
		in[j] = aes_key[j];
		printf("INFO: aes_key[%u] %u\n", j, aes_key[j]);
		// printf("INFO: raw_encrypt_key[%u][%u] | %x\n", indx, j, in[j]);
	}

	for (i = 0; i < val_words; i++, j++)
	{
		in[j] = aes_val[i];
		printf("INFO: aes_val[%u] %u\n", i, aes_val[i]);
		// printf("INFO: raw_encrypt_plaintext[%u][%u] | inputs[%u]@%p %x\n", indx, i, j, in + j, in[j]);
	}

	/* for (j = 0; j < ecb_raw_encrypt_ciphertext_words[indx]; j++) */
	/* { */
	/* 	gold[j] = ecb_raw_encrypt_ciphertext[indx][j]; */
	/* 	// printf("INFO: raw_encrypt_ciphertext[%u][%u] %x\n", indx, j, gold[j]); */
	/* } */
}

static void init_parameters_aes_from_brain(int val_n)
{
	in_words = val_n + (key_length / DATA_BITWIDTH);
	out_words = val_n;

	if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0)
	{
		in_words_adj = in_words;
		out_words_adj = out_words;
	}
	else
	{
		in_words_adj = round_up(in_words, DMA_WORD_PER_BEAT(sizeof(token_t)));
		out_words_adj = round_up(out_words, DMA_WORD_PER_BEAT(sizeof(token_t)));
	}

	in_len = in_words_adj;
	out_len = out_words_adj;
	in_size = in_len * sizeof(token_t);
	out_size = out_len * sizeof(token_t);
	out_offset = in_len;
	size_bytes = (out_offset * sizeof(token_t)) + out_size;
}

int run_brain_bit_only()
{
	printf("==========================================\n");
	printf("==   Start ** brain_bit **              ==\n");
	printf("==========================================\n");

	unsigned errors;

	token_t *gold;
	token_t *buf;

	init_parameters_brain_bit();

	buf = (token_t *)esp_alloc(size);
	cfg_brain_bit_000[0].hw_buf = buf;

	gold = malloc(out_size);

	init_buffer_brain_bit(buf, gold);

	printf("\n====== %s ====== parameters: \n", cfg_brain_bit_000[0].devname);
	printf("  .avg = %f\n", avg);
	printf("  .key_length = %d\n", key_length);
	printf("  .std = %f\n", std);
	printf("  .R = %f\n", R);
	printf("  .L = %d\n", L);
	printf("  .key_batch = %d\n", key_batch);
	printf("  .key_num = %d\n", key_num);

	unsigned avg_u = *avg_ptr;
	unsigned std_u = *std_ptr;
	unsigned R_u = *R_ptr;

	((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->avg = avg_u;
	((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->std = std_u;
	((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->R = R_u;

	// printf("Value avg: %x\n", ((struct brain_bit_vivado_access*) cfg_brain_bit_000[0].esp_desc)->avg);
	// printf("Value std: %x\n", ((struct brain_bit_vivado_access*) cfg_brain_bit_000[0].esp_desc)->std);
	// printf("Value R: %x\n", ((struct brain_bit_vivado_access*) cfg_brain_bit_000[0].esp_desc)->R);

	// token_t* out_location = &buf[out_offset];
	// for(int i = 0; i < out_len; i++)
	// 	out_location[i] = 0;

	printf("\n  ** acc START **\n");
	esp_run(cfg_brain_bit_000, NACC);
	printf("\n  ** acc DONE **\n");

	errors = validate_buffer_brain_bit(&buf[out_offset], gold);

	free(gold);
	esp_free(buf);

	float total = 100 * (float)errors / (key_length * key_batch);
	if (total <= 1)
	{
		printf("- Keys generated: %d\n", key_counter);
		printf("- Number of bit errors: %d\n", errors);
		printf("+ Test PASSED\n\n");
	}
	else
	{
		printf("+ Test FAILED with %f%% errors\n\n", total);
	}

	printf("==========================================\n");
	printf("==   Finish ** brain_bit **             ==\n");
	printf("==========================================\n");

	return errors;
}

int run_aes_only(int n_batch)
{
	printf("==========================================\n");
	printf("==   Start ** aes **                    ==\n");
	printf("==========================================\n");

	unsigned errors = 0;

	token_t *buf_0;
	token_t *gold;


	printf("INFO:   sizeof(token_t) = %lu\n", sizeof(token_t));

	for (unsigned i = 1; i < n_batch; i++)
	{
		printf("\n---------------------------------> aes batch: %d/%d\n", i, n_batch);
		init_parameters_aes(i);
		/* const uint32_t new_in_size = in_size; */
		/* volatile uint32_t new_in_size_reg = *(volatile uint32_t *) &new_in_size; // = in_size - 16; */

		buf_0 = (token_t *)esp_alloc(size_bytes);
		cfg_aes_000[0].hw_buf = buf_0;

		memset(buf_0, 0, size_bytes);

		gold = malloc(out_size);

		init_buffer_aes(buf_0, gold, &buf_0[out_offset], i);

		// printf("\n====== %s ====== parameters: \n", cfg_aes_000[0].devname);
		// printf("  .oper_mode   = %d\n", oper_mode);
		// printf("  .encryption  = %d\n", encryption);
		// printf("  .key_bytes   = %d\n", key_bytes);
		// printf("  .input_bytes = %d\n", input_bytes);
		// printf("  .iv_bytes    = %d\n", iv_bytes);
		// printf("  .aad_bytes   = %d\n", aad_bytes);
		// printf("  .tag_bytes   = %d\n", tag_bytes);
		// printf("  .batch       = %d\n", batch);

		((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->input_bytes = array_size[i]; // new_in_size_reg; //in_size - 16;

		// printf("oper_mode %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->oper_mode);
		// printf("encryption %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->encryption);
		// printf("key_bytes  %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->key_bytes);
		// printf("input_bytes %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->input_bytes);
		// printf("iv_bytes %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->iv_bytes);
		// printf("aad_bytes %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->aad_bytes);
		// printf("tag_bytes %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->tag_bytes);

		printf("\n  ** acc START **\n");
		esp_run(cfg_aes_000, NACC);
		printf("\n  ** acc DONE **\n");

		errors = validate_buffer_aes(&buf_0[out_offset], gold, i);

		free(gold);
		esp_free(buf_0);

		if (!errors)
			printf("  + TEST PASS\n\n");
		else
			printf("  + TEST FAIL\n\n");
	}
	printf("==========================================\n");
	printf("==   Finish ** aes **                   ==\n");
	printf("==========================================\n");

	return errors;
}

int run_brain_bit_and_aes()
{
	printf("==========================================\n");
	printf("==   Start ** brain_bit + aes **        ==\n");
	printf("==========================================\n");

	unsigned errors = 0;

	token_t *gold_brain;
	token_t *buf_brain;
	token_t *buf_aes;

	// set brain_bit parameters
	key_length = 128;
	key_batch = 20;
	key_num = N_TESTS;
	key_num = 1;
	val_num = key_num * 8;

	init_parameters_brain_bit();

	buf_brain = (token_t *)esp_alloc(size);
	cfg_brain_bit_000[0].hw_buf = buf_brain;

	gold_brain = malloc(out_size);

	init_buffer_brain_bit(buf_brain, gold_brain);

	unsigned avg_u = *avg_ptr;
	unsigned std_u = *std_ptr;
	unsigned R_u = *R_ptr;

	((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->avg = avg_u;
	((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->std = std_u;
	((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->R = R_u;

	((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->key_length = key_length;
	((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->key_batch = key_batch;
	((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->key_num = key_num;

	printf("\n====== %s ====== parameters: \n", cfg_brain_bit_000[0].devname);
	printf("  .avg = %f\n", avg);
	printf("  .key_length = %d\n", key_length);
	printf("  .std = %f\n", std);
	printf("  .R = %f\n", R);
	printf("  .L = %d\n", L);
	printf("  .key_batch = %d\n", key_batch);
	printf("  .key_num = %d\n", key_num);

	printf("\n  ** START **\n");
	// Run brain_bit in isolation
	esp_run(cfg_brain_bit_000, NACC);
	printf("\n  ** DONE **\n");

	// set aes keys to key_length_words * key_num
	// token_t aes_key[4*4];
	// set aes plaintext to val_num
	// token_t aes_plain[4*8];

	init_parameters_aes_from_brain(val_num / N_TESTS);

	buf_aes = (token_t *)esp_alloc(size_bytes);
	cfg_aes_000[0].hw_buf = buf_aes;

	memset(buf_aes, 0, size_bytes);

	set_aes_in_from_brain_bit_out(buf_aes, &buf_brain[out_offset]);

	/* init_buffer_aes_from_brain(buf_aes, &buf_brain[out_offset], i); */

	int in_bytes_aes = (val_num / N_TESTS) * sizeof(token_t);

	((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->input_bytes = in_bytes_aes;

	printf("\n====== %s ====== parameters: \n", cfg_aes_000[0].devname);
	printf("  .oper_mode = %d\n", oper_mode);
	printf("  .encryption = %d\n", encryption);
	printf("  .key_bytes = %d\n", key_bytes);
	printf("  .input_bytes = %d\n", input_bytes);
	printf("  .iv_bytes = %d\n", iv_bytes);
	printf("  .aad_bytes = %d\n", aad_bytes);
	printf("  .tag_bytes = %d\n", tag_bytes);
	printf("  .batch = %d\n", batch);
	printf("\n  ** START **\n");

	esp_run(cfg_aes_000, NACC);

	printf("\n  ** DONE **\n");

	// errors = validate_buffer_aes(&buf_aes[out_offset], gold, 0);
	errors = validate_buffer_aes(&buf_aes[out_offset], gold_brain, 0);

	esp_free(buf_brain);

	if (!errors)
		printf("  + TEST PASS\n");
	else
		printf("  + TEST FAIL\n");

	printf("==========================================\n");
	printf("==   Finish ** brain_bit + aes **       ==\n");
	printf("==========================================\n");

	return errors;
}

int main(int argc, char **argv)
{
	int errors_0 = 0;
	int errors_1 = 0;
	int errors_2 = 0;

	errors_0 = run_brain_bit_only();

	errors_1 = run_aes_only(N_BATCH);

	errors_2 = run_brain_bit_and_aes();

	return (errors_0 + errors_1 + errors_2);
}

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

/* User-defined code */
static int validate_buffer_brain_bit(token_t *out, token_t *gold)
{
	int i;
	int j;
	unsigned errors = 0;

	int skip = 0;
	key_counter = 0;
	int val_counter = 0;

	for (i = 0; i < key_batch; i++)
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
						printf("Calculated value %x Golden value %d for index %d \n",
							   val, gold_val, (index - skip));
						if (gold_val != val)
						{
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
					printf("\n----------KEY %d DONE----------\n", key_counter);
					printf("\nKEY IS: [ ");
					for (int k = key_length / 32 - 1; k >= 0; k--)
						printf("0x%x ", out[word - k]);
					printf("]\n\n");
				}
			}
			else if (val_counter != val_num)
			{
				unsigned index = i * out_words_adj + j - skip;
				token_t val = out[index];
				token_t gold_val = (token_t)float_to_fixed32(val_arr[index], 12);
				if (val != gold_val)
					printf("Calculated value %x Golden value %x for index %d \n", val, gold_val, index);
				if ((index + 1) % key_length == 0 && index != 0)
					val_counter++;
			}
		}

	return errors;
}

/* User-defined code */
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

/* User-defined code */
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

/* User-defined code */
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

/* User-defined code */
static void init_buffer_aes(token_t *in, token_t *gold, token_t *out, unsigned indx)
{
	int i;
	int j;

	printf("INFO: raw_encrypt_plaintext_words[%u] %u\n", indx, ecb_raw_encrypt_plaintext_words[indx]);
	printf("INFO: ecb_raw_encrypt_ciphertext_words[%u] %u\n", indx, ecb_raw_encrypt_ciphertext_words[indx]);

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

/* User-defined code */
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
	printf("%s: in_len = %u\n", __func__, in_len);
	out_len = out_words_adj;
	printf("%s: out_len = %u\n", __func__, out_len);
	in_size = in_len * sizeof(token_t);
	printf("%s: in_size = %u\n", __func__, in_size);
	out_size = out_len * sizeof(token_t);
	printf("%s: out_size = %u\n", __func__, out_size);
	out_offset = in_len;
	printf("%s: out_offset = %u\n", __func__, out_offset);
	size_bytes = (out_offset * sizeof(token_t)) + out_size;
	printf("%s: size = %u\n", __func__, size_bytes);
}

int main(int argc, char **argv)
{
	printf("==============================================\n");
	printf("== Start of brain_bit.                      ==\n");
	printf("==============================================\n");

	int errors;

	token_t *gold;
	token_t *buf;

	init_parameters_brain_bit();

	buf = (token_t *)esp_alloc(size);
	cfg_brain_bit_000[0].hw_buf = buf;

	gold = malloc(out_size);

	init_buffer_brain_bit(buf, gold);

	printf("\n====== %s ======\n\n", cfg_brain_bit_000[0].devname);
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

	((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->avg = avg_u;
	((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->std = std_u;
	((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->R = R_u;

	/* printf("Value avg: %x\n", ((struct brain_bit_vivado_access*) cfg_brain_bit_000[0].esp_desc)->avg); */
	/* printf("Value std: %x\n", ((struct brain_bit_vivado_access*) cfg_brain_bit_000[0].esp_desc)->std); */
	/* printf("Value R: %x\n", ((struct brain_bit_vivado_access*) cfg_brain_bit_000[0].esp_desc)->R); */

	/* token_t* out_location = &buf[out_offset]; */
	/* for(int i = 0; i < out_len; i++) */
	/* 	out_location[i] = 0; */

	// esp_run(cfg_brain_bit_000, NACC);

	printf("\n  ** DONE **\n");

	errors = validate_buffer_brain_bit(&buf[out_offset], gold);

	free(gold);
	esp_free(buf);

	float total = 100 * (float)errors / (key_length * key_batch);
	if (total <= 1)
	{
		printf("- Keys generated: %d\n", key_counter);
		printf("- Number of bit errors: %d\n", errors);
		printf("+ Test PASSED\n");
	}
	else
	{
		printf("+ Test FAILED with %f%% errors\n", total);
	}

	printf("\n====== %s ======\n\n", cfg_brain_bit_000[0].devname);

	printf("==============================================\n");
	printf("== End of brain_bit. Start aes.             ==\n");
	printf("==============================================\n");

	unsigned errors_0 = 0;
	unsigned errors_1 = 0;

	// token_t *gold;
	token_t *buf_0;

	printf("INFO:   sizeof(token_t) = %lu\n", sizeof(token_t));

	for (unsigned i = 1; i < N_BATCH; i++)
	{
		init_parameters_aes(i);
		/* const uint32_t new_in_size = in_size; */
		/* volatile uint32_t new_in_size_reg = *(volatile uint32_t *) &new_in_size; // = in_size - 16; */

		buf_0 = (token_t *)esp_alloc(size_bytes);
		cfg_aes_000[0].hw_buf = buf_0;

		memset(buf_0, 0, size_bytes);

		gold = malloc(out_size);

		init_buffer_aes(buf_0, gold, &buf_0[out_offset], i);

		printf("\n====== %s ======\n\n", cfg_aes_000[0].devname);
		/* <<--print-params-->> */
		printf("  .oper_mode = %d\n", oper_mode);
		printf("  .encryption = %d\n", encryption);
		printf("  .key_bytes = %d\n", key_bytes);
		printf("  .input_bytes = %d\n", input_bytes);
		printf("  .iv_bytes = %d\n", iv_bytes);
		printf("  .aad_bytes = %d\n", aad_bytes);
		printf("  .tag_bytes = %d\n", tag_bytes);
		printf("  .batch = %d\n", batch);
		printf("\n  ** START **\n");

		((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->input_bytes = array_size[i]; // new_in_size_reg; //in_size - 16;

		printf("oper_mode %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->oper_mode);
		printf("encryption %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->encryption);
		printf("key_bytes  %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->key_bytes);
		printf("input_bytes %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->input_bytes);
		printf("iv_bytes %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->iv_bytes);
		printf("aad_bytes %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->aad_bytes);
		printf("tag_bytes %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->tag_bytes);

		struct timeval hw_begin_0, hw_end_0;
		gettimeofday(&hw_begin_0, NULL);
		esp_run(cfg_aes_000, NACC);
		gettimeofday(&hw_end_0, NULL);

		printf("\n  ** DONE **\n");

		errors_0 = validate_buffer_aes(&buf_0[out_offset], gold, i);

		free(gold);
		esp_free(buf_0);

		if (!errors_0)
			printf("  + TEST PASS\n");
		else
			printf("  + TEST FAIL\n");

		printf("\n====== %s ======\n\n", cfg_aes_000[0].devname);
	}

	printf("==============================================\n");
	printf("== End of aes.                              ==\n");
	printf("==============================================\n");

	return (errors + errors_0 + errors_1);
}

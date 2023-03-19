// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "libesp.h"
#include "cfg.h"
#include "cfg_p2p_1x1.h"
#include "cfg_p2p_1x2.h"
#include "cfg_p2p_1x3.h"
#include "cfg_p2p_1x4.h"
#include "cfg_p2p_2x1.h"
#include "cfg_p2p_2x3.h"
#include "cfg_p2p_3x1.h"
#include "cfg_p2p_4x1.h"
#include "cfg_p2p_4x4.h"
#include "cfg_p2p_3x2.h"

// #include "brain_bit_input_full.h"
#include "input_1mil_full.h"
#include "aes_data.h"

#include <sys/time.h>

// input_size for each iteration
const unsigned array_size[12] = {16, 32, 48, 64, 80, 96, 112, 128, 144};

void print_brain_bit_config(esp_thread_info_t *x)
{
    printf("\n====== %s ====== config registers: \n", x->devname);
    printf("  .avg          = %d\n", ((struct brain_bit_vivado_access *)x->esp_desc)->avg);
    printf("  .key_length   = %d\n", ((struct brain_bit_vivado_access *)x->esp_desc)->key_length);
    printf("  .std          = %d\n", ((struct brain_bit_vivado_access *)x->esp_desc)->std);
    printf("  .R            = %d\n", ((struct brain_bit_vivado_access *)x->esp_desc)->R);
    printf("  .L            = %d\n", ((struct brain_bit_vivado_access *)x->esp_desc)->L);
    printf("  .key_batch    = %d\n", ((struct brain_bit_vivado_access *)x->esp_desc)->key_batch);
    printf("  .key_num      = %d\n", ((struct brain_bit_vivado_access *)x->esp_desc)->key_num);
    printf("  .val_num      = %d\n", ((struct brain_bit_vivado_access *)x->esp_desc)->val_num);
    printf("  .tot_iter 	= %d\n", ((struct brain_bit_vivado_access *)x->esp_desc)->tot_iter);
}

void print_aes_config(esp_thread_info_t *x)
{
    printf("\n====== %s ====== config registers: \n", x->devname);
    printf("  .oper_mode    = %d\n", ((struct aes_cxx_catapult_access *)x->esp_desc)->oper_mode);
    printf("  .encryption   = %d\n", ((struct aes_cxx_catapult_access *)x->esp_desc)->encryption);
    printf("  .key_bytes    = %d\n", ((struct aes_cxx_catapult_access *)x->esp_desc)->key_bytes);
    printf("  .input_bytes  = %d\n", ((struct aes_cxx_catapult_access *)x->esp_desc)->input_bytes);
    printf("  .iv_bytes     = %d\n", ((struct aes_cxx_catapult_access *)x->esp_desc)->iv_bytes);
    printf("  .aad_bytes    = %d\n", ((struct aes_cxx_catapult_access *)x->esp_desc)->aad_bytes);
    printf("  .tag_bytes    = %d\n", ((struct aes_cxx_catapult_access *)x->esp_desc)->tag_bytes);
    printf("  .batch        = %d\n", ((struct aes_cxx_catapult_access *)x->esp_desc)->batch);
}

void init_parameters_brain_bit(void)
{
    if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0)
    {
        brain_in_words_adj = key_length;
        brain_out_words_adj = key_length;
    }
    else
    {
        brain_in_words_adj = round_up(key_length, DMA_WORD_PER_BEAT(sizeof(token_t)));
        brain_out_words_adj = round_up(key_length, DMA_WORD_PER_BEAT(sizeof(token_t)));
    }
    brain_in_len = brain_in_words_adj * (key_batch);
    brain_out_len = brain_out_words_adj * (key_batch);
    brain_in_size = brain_in_len * sizeof(token_t);
    brain_out_size = brain_out_len * sizeof(token_t);
    brain_out_offset = brain_in_len;
    brain_size = (brain_out_offset * sizeof(token_t)) + brain_out_size;

    printf("%s: brain_in_words_adj  = %u\n", __func__, brain_in_words_adj);
    printf("%s: brain_out_words_adj = %u\n", __func__, brain_out_words_adj);
    printf("%s: brain_in_len        = %u\n", __func__, brain_in_len);
    printf("%s: brain_out_len       = %u\n", __func__, brain_out_len);
    printf("%s: brain_in_size       = %u\n", __func__, brain_in_size);
    printf("%s: brain_out_size      = %u\n", __func__, brain_out_size);
    printf("%s: brain_out_offset    = %u\n", __func__, brain_out_offset);
    printf("%s: brain_size          = %u\n", __func__, brain_size);
}

void init_buffer_brain_bit(token_t *in, token_t *gold)
{
    int i;
    int j;
    int k;

    for (k = 0; k < tot_iter; k++)
    {
        for (i = 0; i < key_batch; i++)
        {
            for (j = 0; j < key_length; j++)
            {
                // printf(" debug 6: k = %d, i = %d, j = %d , index = %d\n", k, i, j, k * key_batch * key_length + i * brain_in_words_adj + j);

                float val = val_arr[k * key_batch * key_length + i * brain_in_words_adj + j];
                // printf(" debug 7: k = %d, i = %d, j = %d , index = %d\n", k, i, j, k * key_batch * key_length + i * brain_in_words_adj + j);

                in[k * key_batch * key_length + i * brain_in_words_adj + j] = (token_t)float_to_fixed32(val, 12);
                // in[i * in_words_adj + j] = (token_t) val;
                /* printf("Generated value %f\n", fixed32_to_float(in[i * out_words_adj + j] , 12)); */
            }
        }
    }

    for (i = 0; i < key_batch; i++)
    {
        for (j = 0; j < key_length; j++)
        {
            // printf(" debug 8: i = %d, j = %d \n", i, j);

            float val = val_arr[i * brain_in_words_adj + j];
            bool filter = (fabs((float)val - avg) >= Rs);
            if (!filter)
            {
                int32_t result = floor((float)(((val - (avg - Rs)) / (2 * Rs)) * L));
                result = result % 2;
                gold[i * brain_out_words_adj + j] = (token_t)result;
                // printf("Generated golden value %d\n", gold[i * out_words_adj + j]);
            }
            else
            {
                gold[i * brain_out_words_adj + j] = 3;
            }
        }
    }
}

int validate_buffer_brain_bit(token_t *out, token_t *gold)
{
    int i;
    int j;
    unsigned errors = 0;

    int skip = 0;
    brain_key_counter = 0;
    /* int val_counter = 0; */
    int offset = 0;
    bool done = false;

    for (i = 0; i < key_batch; i++)
    {
        for (j = 0; j < key_length; j++)
        {
            if (brain_key_counter != key_num)
            {
                unsigned index = i * brain_out_words_adj + j;
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

                if ((index - skip + 1) % (key_length * (brain_key_counter + 1)) == 0 && index != 0)
                {
                    brain_key_counter++;
                    printf("----------KEY %d DONE----------\n", brain_key_counter);
                    printf("KEY IS: [ ");
                    for (int k = key_length / 32 - 1; k >= 0; k--)
                        printf("0x%x ", out[word - k]);
                    printf("]\n");
                }
            }
            else if (!done)
            {
                done = true;
                offset = i * brain_out_words_adj + j - skip;
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

void init_parameters_aes(unsigned indx)
{
    aes_in_words = ecb_plaintext_words[indx] +
                   ecb_key_words[indx];
    aes_out_words = ecb_plaintext_words[indx];

    if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0)
    {
        aes_in_words_adj = aes_in_words;
        aes_out_words_adj = aes_out_words;
    }
    else
    {
        aes_in_words_adj = round_up(aes_in_words, DMA_WORD_PER_BEAT(sizeof(token_t)));
        aes_out_words_adj = round_up(aes_out_words, DMA_WORD_PER_BEAT(sizeof(token_t)));
    }

    aes_in_len = aes_in_words_adj;
    aes_out_len = aes_out_words_adj;
    aes_in_size = aes_in_len * sizeof(token_t);
    aes_out_size = aes_out_len * sizeof(token_t);
    aes_out_offset = aes_in_len;
    aes_size_bytes = (aes_out_offset * sizeof(token_t)) + aes_out_size;

    printf("%s: aes_in_words_adj  = %u\n", __func__, aes_in_words_adj);
    printf("%s: aes_out_words_adj = %u\n", __func__, aes_out_words_adj);
    printf("%s: aes_in_len        = %u\n", __func__, aes_in_len);
    printf("%s: aes_out_len       = %u\n", __func__, aes_out_len);
    printf("%s: aes_in_size       = %u\n", __func__, aes_in_size);
    printf("%s: aes_out_size      = %u\n", __func__, aes_out_size);
    printf("%s: aes_out_offset    = %u\n", __func__, aes_out_offset);
    printf("%s: aes_size_bytes    = %u\n", __func__, aes_size_bytes);
}

void init_buffer_aes(token_t *in, token_t *gold, token_t *out, unsigned indx)
{
    int i;
    int j;

    printf("init_buffer_aes: ecb_plaintext_words [%u] %u\n", indx, ecb_plaintext_words[indx]);
    printf("init_buffer_aes: ecb_ciphertext_words[%u] %u\n", indx, ecb_ciphertext_words[indx]);

    for (j = 0; j < ecb_key_words[indx]; j++)
    {
        in[j] = ecb_key[indx][j];
        // printf("INFO: raw_encrypt_key[%u][%u] | %x\n", indx, j, in[j]);
    }

    for (i = 0; i < ecb_plaintext_words[indx]; i++, j++)
    {
        in[j] = ecb_plaintext[indx][i];
        // printf("INFO: raw_encrypt_plaintext[%u][%u] | inputs[%u]@%p %x\n", indx, i, j, in + j, in[j]);
    }

    for (j = 0; j < ecb_ciphertext_words[indx]; j++)
    {
        gold[j] = ecb_ciphertext[indx][j];
        // printf("INFO: raw_encrypt_ciphertext[%u][%u] %x\n", indx, j, gold[j]);
    }
}

void init_buffer_aes_p2p_1x1(token_t *in, token_t *gold, unsigned indx)
{
    // int i;
    int j;

    printf("init_buffer_aes_p2p_1x1: ecb_plaintext_words_p2p_1x1  [%u] %u\n", indx, ecb_plaintext_words_p2p_1x1[indx]);
    printf("init_buffer_aes_p2p_1x1: ecb_ciphertext_words_p2p_1x1 [%u] %u\n", indx, ecb_ciphertext_words_p2p_1x1[indx]);

    // for (j = 0; j < ecb_key_words_p2p_1x1[indx]; j++)
    // {
    //     in[j] = ecb_key_p2p_1x1[indx][j];
    //     printf("INFO: raw_encrypt_key_p2p_1x1[%u][%u] | %x\n", indx, j, in[j]);
    // }

    // for (i = 0; i < ecb_plaintext_words_p2p_1x1[indx]; i++, j++)
    // {
    //     in[j] = ecb_plaintext_p2p_1x1[indx][i];
    //     printf("INFO: raw_encrypt_plaintext_p2p_1x1[%u][%u] | inputs[%u]@%p %x\n", indx, i, j, in + j, in[j]);
    // }

    for (j = 0; j < ecb_ciphertext_words_p2p_1x1[indx]; j++)
    {
        gold[j] = ecb_ciphertext_p2p_1x1[indx][j];
        printf("INFO: raw_encrypt_ciphertext_p2p_1x1[%u][%u] %x\n", indx, j, gold[j]);
    }
}

int validate_buffer_aes(token_t *in, token_t *out, token_t *gold, unsigned indx)
{
    /* int i; */
    int j;
    unsigned errors = 0;

    printf("INFO:  gold output data @%p\n", gold);
    printf("INFO:       output data @%p\n", out);

    for (j = 0; j < ecb_ciphertext_words[indx] /* + 0x10 */; j++)
    {
        token_t in_data = in[j];
        token_t gold_data = gold[j];
        token_t out_data = out[j];

        if (out_data != gold_data)
        {
            errors++;
            printf("mistmatch_aes: [%d] @%p, in: %x, out: %x, gold: %x\n", j, out + j, in_data, out_data, gold_data);
        }
        // printf("validate_aes: [%d] @%p, in: %x, out: %x, gold: %x\n", j, out + j, in_data, out_data, gold_data);
    }

    return errors;
}

void set_aes_in_from_brain_bit_out(token_t *in_aes, token_t *out_brain)
{
    // this function should only take 1 key

    int i;
    int j = 0;

    int key_num_words = key_num * key_length / DATA_BITWIDTH;
    /* int val_num_words = val_num * key_length / DATA_BITWIDTH; */
    int in_num_words = (key_length - 128) / DATA_BITWIDTH;

    if (DMA_WORD_PER_BEAT(sizeof(token_t)) != 0)
    {
        key_num_words = round_up(key_num_words, DMA_WORD_PER_BEAT(sizeof(token_t)));
        in_num_words = round_up(in_num_words, DMA_WORD_PER_BEAT(sizeof(token_t)));
    }

    for (i = 0; i < key_num_words; i++)
    {
        in_aes[j] = out_brain[i];
	if(i > 0 && i % (key_length / DATA_BITWIDTH - 1) == 0)
		j += in_num_words;//jump to next input
	else
		j++;

        // printf("in_aes[%d] = %x\n", i, in_aes[i]);
    }
    printf("---\n");
    /* for (j = 0; j < val_num_words; j++) */
    /* { */
    /*     in_aes[i + j] = out_brain[j + i]; */
    /*     // printf("in_aes[%d] = %x\n", i + j, in_aes[i + j]); */
    /* } */
}

/*
void init_buffer_aes_from_brain(token_t *in, token_t *aes_key, token_t *aes_val, token_t *out, unsigned indx)
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

    for (j = 0; j < ecb_ciphertext_words[indx]; j++)
    {
        gold[j] = ecb_ciphertext[indx][j];
        // printf("INFO: raw_encrypt_ciphertext[%u][%u] %x\n", indx, j, gold[j]);
    }
}
*/

void init_parameters_aes_from_brain()
{
    aes_in_words = key_length*key_num / DATA_BITWIDTH;
    aes_out_words = (key_length*key_num-128) / DATA_BITWIDTH;

    if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0)
    {
        aes_in_words_adj = aes_in_words;
        aes_out_words_adj = aes_out_words;
    }
    else
    {
        aes_in_words_adj = round_up(aes_in_words, DMA_WORD_PER_BEAT(sizeof(token_t)));
        aes_out_words_adj = round_up(aes_out_words, DMA_WORD_PER_BEAT(sizeof(token_t)));
    }

    aes_in_len = aes_in_words_adj;
    aes_out_len = aes_out_words_adj;
    aes_in_size = aes_in_len * sizeof(token_t);
    aes_out_size = aes_out_len * sizeof(token_t);
    aes_out_offset = aes_in_len;
    aes_size_bytes = (aes_out_offset * sizeof(token_t)) + aes_out_size;
}

int run_brain_bit_only(void)
{
    printf("==================================================\n");
    printf("==   Start ** brain_bit **                      ==\n");
    printf("==================================================\n");

    unsigned errors;

    token_t *gold;
    token_t *buf;

    init_parameters_brain_bit();

    buf = (token_t *)esp_alloc(brain_size);
    cfg_brain_bit_000[0].hw_buf = buf;
    memset(buf, 0, brain_size);

    gold = malloc(brain_out_size);

    init_buffer_brain_bit(buf, gold);

    printf("\n====== %s ====== config registers: \n", cfg_brain_bit_000[0].devname);
    printf("  .avg          = %f\n", avg);
    printf("  .key_length   = %d\n", key_length);
    printf("  .std          = %f\n", std);
    printf("  .R            = %f\n", R);
    printf("  .L            = %d\n", L);
    printf("  .key_batch    = %d\n", key_batch);
    printf("  .key_num      = %d\n", key_num);
    printf("  .val_num      = %d\n", val_num);
    printf("  .tot_iter 	= %d\n", tot_iter);

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

    errors = validate_buffer_brain_bit(&buf[brain_out_offset], gold);

    float total = 100 * (float)errors / (key_length * key_batch);
    if (total <= 1)
    {
        printf("- Keys generated: 		%d\n", brain_key_counter);
        printf("- Number of bit errors: %d\n", errors);
        printf("+ Test PASSED\n\n");
    }
    else
    {
        printf("+ Test FAILED with %f%% errors\n\n", total);
    }

    printf("==================================================\n");
    printf("==   Finish ** brain_bit **                     ==\n");
    printf("==================================================\n");

    free(gold);
    esp_free(buf);

    return errors;
}

int run_aes_only(int n_batch)
{
    printf("==================================================\n");
    printf("==   Start ** aes **                            ==\n");
    printf("==================================================\n");

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

        buf_0 = (token_t *)esp_alloc(aes_size_bytes);
        cfg_aes_000[0].hw_buf = buf_0;
        memset(buf_0, 0, aes_size_bytes);

        gold = malloc(aes_out_size);

        init_buffer_aes(buf_0, gold, &buf_0[aes_out_offset], i);

        printf("\n====== %s ====== config registers: \n", cfg_aes_000[0].devname);
        printf("  .oper_mode   = %d\n", oper_mode);
        printf("  .encryption  = %d\n", encryption);
        printf("  .key_bytes   = %d\n", key_bytes);
        printf("  .input_bytes = %d\n", input_bytes);
        printf("  .iv_bytes    = %d\n", iv_bytes);
        printf("  .aad_bytes   = %d\n", aad_bytes);
        printf("  .tag_bytes   = %d\n", tag_bytes);
        printf("  .batch       = %d\n", batch);

        ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->input_bytes = array_size[i]; // new_in_size_reg; //in_size - 16;

        printf("oper_mode   = %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->oper_mode);
        printf("encryption  = %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->encryption);
        printf("key_bytes   = %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->key_bytes);
        printf("input_bytes = %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->input_bytes);
        printf("iv_bytes    = %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->iv_bytes);
        printf("aad_bytes   = %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->aad_bytes);
        printf("tag_bytes   = %u\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->tag_bytes);

        printf("\n  ** acc START **\n");
        esp_run(cfg_aes_000, NACC);
        printf("\n  ** acc DONE **\n");

        errors = validate_buffer_aes(buf_0, &buf_0[aes_out_offset], gold, i);
        printf("INFO: total errors %u\n", errors);

        if (!errors)
            printf("  + TEST PASS\n\n");
        else
            printf("  + TEST FAIL\n\n");

        free(gold);
        esp_free(buf_0);
    }
    printf("==================================================\n");
    printf("==   Finish ** aes **                           ==\n");
    printf("==================================================\n");

    return errors;
}

int run_both_mem_1x1(int N)
{
    printf("==================================================\n");
    printf("==   Start ** brain_bit + aes (mem_1x1) **      ==\n");
    printf("==================================================\n");

    unsigned long long total_time = 0;

    unsigned errors = 0;
    /* int i = 0; */

    token_t *buf_brain;
    /* token_t *gold_brain; */

    token_t *buf_aes;
    /* token_t *gold_aes; */

    // set brain_bit config registers
    //key_length = 128;
    key_batch = N*1.1;
    key_num = N;
    val_num = 0;
    tot_iter = 1;

    init_parameters_brain_bit();
    init_parameters_aes_from_brain();

    buf_brain = (token_t *)esp_alloc(brain_size);
    /* buf_aes = (token_t *)esp_alloc(aes_size_bytes); */
    cfg_brain_bit_000[0].hw_buf = buf_brain;
    memset(buf_brain, 0, brain_size);
    /* gold_brain = malloc(brain_out_size); */

    /* init_buffer_brain_bit(buf_brain, gold_brain); */

    /* gold_brain = malloc(brain_out_size); */
    /* gold_aes = malloc(aes_out_size); */

    unsigned avg_u = *avg_ptr;
    unsigned std_u = *std_ptr;
    unsigned R_u = *R_ptr;

    /* print_brain_bit_config(&cfg_brain_bit_000[0]); */

    ((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_brain_bit_000[0].esp_desc)->tot_iter = tot_iter;

    token_t* brain_out_location = &buf_brain[brain_out_offset];
    /* buf_aes = brain_out_location; */

    print_brain_bit_config(&cfg_brain_bit_000[0]);

    /* printf("\n====== %s ====== config registers: \n", cfg_brain_bit_000[0].devname); */
    /* printf("  .avg          = %f\n", avg); */
    /* printf("  .key_length   = %d\n", key_length); */
    /* printf("  .std          = %f\n", std); */
    /* printf("  .R            = %f\n", R); */
    /* printf("  .L            = %d\n", L); */
    /* printf("  .key_batch    = %d\n", key_batch); */
    /* printf("  .key_num      = %d\n", key_num); */
    /* printf("  .val_num      = %d\n", val_num); */
    /* printf("  .tot_iter 	= %d\n", tot_iter); */

    /* printf("\n  ** START **\n"); */
    /* // Run brain_bit in isolation */
    /* for (i = 0; i < N; i++) */
    /* { */
    /*     esp_run(cfg_brain_bit_000, NACC); */
    /*     total_time += cfg_brain_bit_000[0].hw_ns; */
    /* } */
    /* printf("\n  ** DONE **\n"); */

    // set aes keys to key_length_words * key_num
    // token_t aes_key[4*4];
    // set aes plaintext to val_num
    // token_t aes_plain[4*8];

    /* init_parameters_aes_from_brain(); */
    // init_parameters_aes(1);

    buf_aes = (token_t *)esp_alloc(aes_size_bytes);
    cfg_aes_000[0].hw_buf = buf_aes;
    memset(buf_aes, 0, aes_size_bytes);

    //gold_aes = malloc(aes_out_size);

    /* set_aes_in_from_brain_bit_out(buf_aes, &buf_brain[brain_out_offset]); */
    /* init_buffer_aes_p2p_1x1(buf_aes, gold_aes, 0); */

    // init_buffer_aes_from_brain(buf_aes, &buf_brain[brain_out_offset], i);

    /* input_bytes = ((key_length-128)/DATA_BITWIDTH) * sizeof(token_t); */

    ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->batch = N;

    print_aes_config(&cfg_aes_000[0]);

    /* printf("\n====== %s ====== config registers: \n", cfg_aes_000[0].devname); */
    /* printf("  .oper_mode    = %d\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->oper_mode); */
    /* printf("  .encryption   = %d\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->encryption); */
    /* printf("  .key_bytes    = %d\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->key_bytes); */
    /* printf("  .input_bytes  = %d\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->input_bytes); */
    /* printf("  .iv_bytes     = %d\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->iv_bytes); */
    /* printf("  .aad_bytes    = %d\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->aad_bytes); */
    /* printf("  .tag_bytes    = %d\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->tag_bytes); */
    /* printf("  .batch        = %d\n", ((struct aes_cxx_catapult_access *)cfg_aes_000[0].esp_desc)->batch); */

    struct timespec th_start;
    struct timespec th_end;
    unsigned long long brain_time = 0;
    unsigned long long aes_time = 0;
    unsigned long long mem_time = 0;

    printf("\n  ** START **\n");
    /* for (i = 0; i < N; i++) */
    /* { */

    gettime(&th_start);
    esp_run_no_print(cfg_brain_bit_000, 1);
    gettime(&th_end);

    brain_time = ts_subtract(&th_start, &th_end);

    gettime(&th_start);
    set_aes_in_from_brain_bit_out(buf_aes, brain_out_location);
    gettime(&th_end);

    mem_time = ts_subtract(&th_start, &th_end);

    gettime(&th_start);
    esp_run_no_print(cfg_aes_000, 1);
    gettime(&th_end);

    aes_time = ts_subtract(&th_start, &th_end);

    total_time = cfg_aes_000[0].hw_ns + cfg_brain_bit_000[0].hw_ns;
    /* } */
    printf("\n  ** DONE **\n");

    // errors = validate_buffer_aes(buf_aes, &buf_aes[aes_out_offset], gold_aes, 0);
    printf("INFO: total errors %u\n", errors);

    if (!errors)
        printf("  + TEST PASS\n");
    else
        printf("  + TEST FAIL\n");

    total_time = brain_time+mem_time+aes_time;

    printf("==================================================\n");
    printf("==   Finish ** brain_bit + aes (mem_1x1) **     ==\n");
    printf("==     N: %d\tTotal time: %llu\n", N, total_time);
    printf("==================================================\n");

    fprintf(log_file, "-- mem_1x1 brain_bit time: %llu\n", cfg_brain_bit_000[0].hw_ns);
    fprintf(log_file, "-- mem_1x1 aes time: %llu\n", cfg_aes_000[0].hw_ns);
    fprintf(log_file, "-- mem_1x1 N: %d time:\t%llu\n", N, total_time);

    esp_free(buf_brain);
    esp_free(buf_aes);

    return total_time;
}

int run_both_p2p_1x1(int N)
{
    printf("==================================================\n");
    printf("==   Start ** brain_bit + aes (p2p_1x1) **      ==\n");
    printf("==================================================\n");

    unsigned long long total_time = 0;

    unsigned errors = 0;

    token_t *buf_brain;
    /* token_t *gold_brain; */

    token_t *buf_aes;
    /* token_t *gold_aes; */

    // set brain_bit config registers
    // key_length = 512; // 128;
    key_batch = N*1.1;  // 3;
    key_num = N; // 1;
    val_num = 0; // key_num * 4;
    tot_iter = 1; // N;

    init_parameters_brain_bit();
    // init_parameters_aes_from_brain(val_num);

    // buf_brain = (token_t *)esp_alloc(brain_size * tot_iter);
    // // buf_aes = (token_t *)esp_alloc(aes_size_bytes);
    // buf_aes = (token_t *)esp_alloc(brain_size * tot_iter);
    // cfg_p2p_1x1[0].hw_buf = buf_brain;
    // cfg_p2p_1x1[1].hw_buf = buf_aes;
    // memset(buf_brain, 0, brain_size * tot_iter);
    // // memset(buf_aes, 0, aes_size_bytes);
    // memset(buf_aes, 0, brain_size * tot_iter);

    buf_brain = (token_t *)esp_alloc(brain_size);
    buf_aes = (token_t *)esp_alloc(brain_size);
    cfg_p2p_1x1[0].hw_buf = buf_brain;
    cfg_p2p_1x1[1].hw_buf = buf_aes;
    memset(buf_brain, 0, brain_size);
    memset(buf_aes, 0, brain_size);


    /* gold_brain = malloc(brain_out_size); */
    /* gold_aes = malloc(aes_out_size); */

    // init_buffer_brain_bit(buf_brain, gold_brain);
    unsigned avg_u = *avg_ptr;
    unsigned std_u = *std_ptr;
    unsigned R_u = *R_ptr;

    ((struct brain_bit_vivado_access *)cfg_p2p_1x1[0].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_1x1[0].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_1x1[0].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_1x1[0].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_1x1[0].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_1x1[0].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_1x1[0].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_1x1[0].esp_desc)->tot_iter = tot_iter;

    set_aes_in_from_brain_bit_out(buf_aes, &buf_brain[brain_out_offset]);
    /* init_buffer_aes_p2p_1x1(buf_aes, gold_aes, 0); */

    // input_bytes = 48; // val_num * sizeof(token_t);

    ((struct aes_cxx_catapult_access *)cfg_p2p_1x1[1].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_1x1[1].esp_desc)->batch = N;

    print_brain_bit_config(&cfg_p2p_1x1[0]);
    print_aes_config(&cfg_p2p_1x1[1]);

    printf("\n  ** START **\n");
    esp_run(cfg_p2p_1x1, 2);
    total_time = max(cfg_p2p_1x1[0].hw_ns, cfg_p2p_1x1[1].hw_ns);
    printf("\n  ** DONE **\n");

    // errors = validate_buffer_aes(buf_aes, &buf_aes[aes_out_offset], gold_aes, 0);
    printf("INFO: total errors %u\n", errors);

    if (!errors)
        printf("  + TEST PASS\n");
    else
        printf("  + TEST FAIL\n");

    printf("==================================================\n");
    printf("==   Finish ** brain_bit + aes (p2p_1x1) **     ==\n");
    printf("==     N: %d\tTotal time: %llu\n", N, total_time);
    printf("==================================================\n");

    fprintf(log_file, "-- p2p_1x1[0] time: %llu\n", cfg_p2p_1x1[0].hw_ns);
    fprintf(log_file, "-- p2p_1x1[1] time: %llu\n", cfg_p2p_1x1[1].hw_ns);
    fprintf(log_file, "-- p2p_1x1 N: %d time:\t%llu\n", N, total_time);

    esp_free(buf_brain);
    esp_free(buf_aes);

    return total_time;
}

int run_both_p2p_1x2(int N)
{
    printf("==================================================\n");
    printf("==   Start ** brain_bit + aes (p2p_1x2) **      ==\n");
    printf("==================================================\n");

    unsigned long long total_time = 0;

    unsigned errors = 0;

    token_t *buf_brain;
    /* token_t *gold_brain; */

    token_t *buf_aes;
    /* token_t *gold_aes; */

    // set brain_bit config registers
    // key_length = 128;
    key_batch = 2*N*1.1; // 3;
    key_num = 2*N; // 1
    val_num = 0; // key_num * 4;
    tot_iter = 1; // 2 * N;

    init_parameters_brain_bit();
    // init_parameters_aes_from_brain(val_num);


    buf_brain = (token_t *)esp_alloc(brain_size);
    buf_aes = (token_t *)esp_alloc(brain_size);
    cfg_p2p_1x2[0].hw_buf = buf_brain;
    cfg_p2p_1x2[1].hw_buf = buf_aes;
    cfg_p2p_1x2[2].hw_buf = buf_aes;
    memset(buf_brain, 0, brain_size);
    memset(buf_aes, 0, brain_size);


    /* gold_brain = malloc(brain_out_size); */
    /* gold_aes = malloc(aes_out_size); */

    // init_buffer_brain_bit(buf_brain, gold_brain);
    unsigned avg_u = *avg_ptr;
    unsigned std_u = *std_ptr;
    unsigned R_u = *R_ptr;

    ((struct brain_bit_vivado_access *)cfg_p2p_1x2[0].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_1x2[0].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_1x2[0].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_1x2[0].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_1x2[0].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_1x2[0].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_1x2[0].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_1x2[0].esp_desc)->tot_iter = tot_iter;

    set_aes_in_from_brain_bit_out(buf_aes, &buf_brain[brain_out_offset]);
    /* init_buffer_aes_p2p_1x1(buf_aes, gold_aes, 0); */

    // input_bytes = val_num * sizeof(token_t);

    ((struct aes_cxx_catapult_access *)cfg_p2p_1x2[1].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_1x2[1].esp_desc)->batch = N;
    ((struct aes_cxx_catapult_access *)cfg_p2p_1x2[2].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_1x2[2].esp_desc)->batch = N;

    print_brain_bit_config(&cfg_p2p_1x2[0]);
    print_aes_config(&cfg_p2p_1x2[1]);
    print_aes_config(&cfg_p2p_1x2[2]);

    printf("\n  ** START **\n");
    esp_run(cfg_p2p_1x2, 3);
    total_time = max(max(cfg_p2p_1x2[0].hw_ns, cfg_p2p_1x2[1].hw_ns), cfg_p2p_1x2[2].hw_ns);
    printf("\n  ** DONE **\n");

    // errors = validate_buffer_aes(buf_aes, &buf_aes[aes_out_offset], gold_aes, 0);
    printf("INFO: total errors %u\n", errors);

    if (!errors)
        printf("  + TEST PASS\n");
    else
        printf("  + TEST FAIL\n");

    printf("==================================================\n");
    printf("==   Finish ** brain_bit + aes (p2p_1x2) **     ==\n");
    printf("==================================================\n");

    fprintf(log_file, "-- p2p_1x2[0] time: %llu\n", cfg_p2p_1x2[0].hw_ns);
    fprintf(log_file, "-- p2p_1x2[1] time: %llu\n", cfg_p2p_1x2[1].hw_ns);
    fprintf(log_file, "-- p2p_1x2[2] time: %llu\n", cfg_p2p_1x2[2].hw_ns);
    fprintf(log_file, "-- p2p_1x2 N: %d time:\t%llu\n", N, total_time);

    esp_free(buf_brain);
    esp_free(buf_aes);

    return total_time;
}

int run_both_p2p_2x1(int N)
{
    printf("==================================================\n");
    printf("==   Start ** brain_bit + aes (p2p_2x1) **      ==\n");
    printf("==================================================\n");

    unsigned long long total_time = 0;

    unsigned errors = 0;

    token_t *buf_brain;
    /* token_t *gold_brain; */

    token_t *buf_aes;
    /* token_t *gold_aes; */

    // set brain_bit config registers
    // key_length = 128;
    key_batch = N*1.1;  // 3;
    key_num = N; // 1;
    val_num = 0; // key_num * 4;
    tot_iter = 1; // N;

    init_parameters_brain_bit();
    // init_parameters_aes_from_brain(val_num);

    buf_brain = (token_t *)esp_alloc(brain_size);
    buf_aes = (token_t *)esp_alloc(brain_size);
    cfg_p2p_2x1[0].hw_buf = buf_brain;
    cfg_p2p_2x1[1].hw_buf = buf_brain;
    cfg_p2p_2x1[2].hw_buf = buf_aes;

    printf("aes_size = %d\n", aes_size_bytes);
    printf("aes_out_size = %d\n", aes_out_size);
    printf("brain_size = %d\n", brain_size);
    printf("brain_out_size = %d\n", brain_out_size);

    memset(buf_brain, 0, brain_size);
    memset(buf_aes, 0, brain_size);

    /* gold_brain = malloc(brain_out_size * tot_iter); */
    /* gold_aes = malloc(aes_out_size * tot_iter); */

    // init_buffer_brain_bit(buf_brain, gold_brain);
    unsigned avg_u = *avg_ptr;
    unsigned std_u = *std_ptr;
    unsigned R_u = *R_ptr;

    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[0].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[0].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[0].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[0].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[0].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[0].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[0].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[0].esp_desc)->tot_iter = tot_iter;

    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[1].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[1].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[1].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[1].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[1].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[1].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[1].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x1[1].esp_desc)->tot_iter = tot_iter;

    set_aes_in_from_brain_bit_out(buf_aes, &buf_brain[brain_out_offset]);
    /* init_buffer_aes_p2p_1x1(buf_aes, gold_aes, 0); */

    // input_bytes = val_num * sizeof(token_t);

    ((struct aes_cxx_catapult_access *)cfg_p2p_2x1[2].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_2x1[2].esp_desc)->batch = 2 * N;

    print_brain_bit_config(&cfg_p2p_1x1[0]);
    print_brain_bit_config(&cfg_p2p_1x1[1]);
    print_aes_config(&cfg_p2p_1x1[2]);

    printf("\n  ** START **\n");
    esp_run(cfg_p2p_2x1, 3);
    total_time = max(max(cfg_p2p_2x1[0].hw_ns, cfg_p2p_2x1[1].hw_ns), cfg_p2p_2x1[2].hw_ns);
    printf("\n  ** DONE **\n");

    // errors = validate_buffer_aes(buf_aes, &buf_aes[aes_out_offset], gold_aes, 0);
    printf("INFO: total errors %u\n", errors);

    if (!errors)
        printf("  + TEST PASS\n");
    else
        printf("  + TEST FAIL\n");

    printf("==================================================\n");
    printf("==   Finish ** brain_bit + aes (p2p_2x1) **     ==\n");
    printf("==     N: %d\tTotal time: %llu\n", N, total_time);
    printf("==================================================\n");

    fprintf(log_file, "-- p2p_2x1[0] time: %llu\n", cfg_p2p_2x1[0].hw_ns);
    fprintf(log_file, "-- p2p_2x1[1] time: %llu\n", cfg_p2p_2x1[1].hw_ns);
    fprintf(log_file, "-- p2p_2x1[2] time: %llu\n", cfg_p2p_2x1[2].hw_ns);
    fprintf(log_file, "-- p2p_2x1 N: %d time:\t%llu\n", N, total_time);

    esp_free(buf_brain);
    esp_free(buf_aes);

    return total_time;
}

int run_both_p2p_3x1(int N)
{
    printf("==================================================\n");
    printf("==   Start ** brain_bit + aes (p2p_3x1) **      ==\n");
    printf("==================================================\n");

    unsigned long long total_time = 0;

    unsigned errors = 0;

    token_t *buf_brain;
    /* token_t *gold_brain; */

    token_t *buf_aes;
    /* token_t *gold_aes; */

    // set brain_bit config registers
    // key_length = 128;
    // key_batch = 3;
    // key_num = 1;
    // val_num = 0; // key_num * 4;
    // tot_iter = N;
    key_batch = N*1.1;
    key_num = N;
    val_num = 0;
    tot_iter = 1;

    init_parameters_brain_bit();
    // init_parameters_aes_from_brain(val_num);

    buf_brain = (token_t *)esp_alloc(brain_size);
    buf_aes = (token_t *)esp_alloc(brain_size);
    cfg_p2p_3x1[0].hw_buf = buf_brain;
    cfg_p2p_3x1[1].hw_buf = buf_brain;
    cfg_p2p_3x1[2].hw_buf = buf_brain;
    cfg_p2p_3x1[3].hw_buf = buf_aes;
    memset(buf_brain, 0, brain_size);
    memset(buf_aes, 0, brain_size);

    /* gold_brain = malloc(brain_out_size * tot_iter); */
    /* gold_aes = malloc(aes_out_size * tot_iter); */

    // init_buffer_brain_bit(buf_brain, gold_brain);

    unsigned avg_u = *avg_ptr;
    unsigned std_u = *std_ptr;
    unsigned R_u = *R_ptr;

    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[0].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[0].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[0].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[0].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[0].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[0].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[0].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[0].esp_desc)->tot_iter = tot_iter;

    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[1].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[1].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[1].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[1].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[1].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[1].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[1].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[1].esp_desc)->tot_iter = tot_iter;

    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[2].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[2].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[2].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[2].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[2].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[2].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[2].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x1[2].esp_desc)->tot_iter = tot_iter;

    set_aes_in_from_brain_bit_out(buf_aes, &buf_brain[brain_out_offset]);
    /* init_buffer_aes_p2p_1x1(buf_aes, gold_aes, 0); */

    // input_bytes = val_num * sizeof(token_t);

    ((struct aes_cxx_catapult_access *)cfg_p2p_3x1[3].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_3x1[3].esp_desc)->batch = 3 * N;

    print_brain_bit_config(&cfg_p2p_3x1[0]);
    print_brain_bit_config(&cfg_p2p_3x1[1]);
    print_brain_bit_config(&cfg_p2p_3x1[2]);
    print_aes_config(&cfg_p2p_3x1[3]);

    printf("\n  ** START **\n");
    esp_run(cfg_p2p_3x1, 4);
    total_time = max(max(cfg_p2p_3x1[0].hw_ns,
                         max(cfg_p2p_3x1[1].hw_ns,
                             cfg_p2p_3x1[2].hw_ns)),
                     cfg_p2p_3x1[3].hw_ns);
    printf("\n  ** DONE **\n");

    // errors = validate_buffer_aes(buf_aes, &buf_aes[aes_out_offset], gold_aes, 0);
    printf("INFO: total errors %u\n", errors);

    if (!errors)
        printf("  + TEST PASS\n");
    else
        printf("  + TEST FAIL\n");

    printf("==================================================\n");
    printf("==   Finish ** brain_bit + aes (p2p_3x1) **     ==\n");
    printf("==     N: %d\tTotal time: %llu\n", N, total_time);
    printf("==================================================\n");

    fprintf(log_file, "-- p2p_3x1[0] time: %llu\n", cfg_p2p_3x1[0].hw_ns);
    fprintf(log_file, "-- p2p_3x1[1] time: %llu\n", cfg_p2p_3x1[1].hw_ns);
    fprintf(log_file, "-- p2p_3x1[2] time: %llu\n", cfg_p2p_3x1[2].hw_ns);
    fprintf(log_file, "-- p2p_3x1[3] time: %llu\n", cfg_p2p_3x1[3].hw_ns);
    fprintf(log_file, "-- p2p_3x1 N: %d time:\t%llu\n", N, total_time);

    esp_free(buf_brain);
    esp_free(buf_aes);

    return total_time;
}

int run_both_p2p_4x1(int N)
{
    printf("==================================================\n");
    printf("==   Start ** brain_bit + aes (p2p_4x1) **      ==\n");
    printf("==================================================\n");

    unsigned long long total_time = 0;

    unsigned errors = 0;

    token_t *buf_brain;
    /* token_t *gold_brain; */

    token_t *buf_aes;
    /* token_t *gold_aes; */

    // set brain_bit config registers
    // key_length = 128;
    // key_batch = 3;
    // key_num = 1;
    // val_num = 0; // key_num * 4;
    // tot_iter = N;
    key_batch = N*1.1;
    key_num = N;
    val_num = 0;
    tot_iter = 1;

    init_parameters_brain_bit();
    // init_parameters_aes_from_brain(val_num);

    buf_brain = (token_t *)esp_alloc(brain_size);
    buf_aes = (token_t *)esp_alloc(brain_size);
    cfg_p2p_4x1[0].hw_buf = buf_brain;
    cfg_p2p_4x1[1].hw_buf = buf_brain;
    cfg_p2p_4x1[2].hw_buf = buf_brain;
    cfg_p2p_4x1[3].hw_buf = buf_brain;
    cfg_p2p_4x1[4].hw_buf = buf_aes;
    memset(buf_brain, 0, brain_size);
    memset(buf_aes, 0, brain_size);

    /* gold_brain = malloc(brain_out_size * N); */
    /* gold_aes = malloc(aes_out_size); */

    // init_buffer_brain_bit(buf_brain, gold_brain);

    unsigned avg_u = *avg_ptr;
    unsigned std_u = *std_ptr;
    unsigned R_u = *R_ptr;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[0].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[0].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[0].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[0].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[0].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[0].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[0].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[0].esp_desc)->tot_iter = tot_iter;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[1].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[1].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[1].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[1].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[1].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[1].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[1].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[1].esp_desc)->tot_iter = tot_iter;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[2].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[2].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[2].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[2].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[2].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[2].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[2].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[2].esp_desc)->tot_iter = tot_iter;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[3].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[3].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[3].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[3].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[3].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[3].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[3].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x1[3].esp_desc)->tot_iter = tot_iter;

    set_aes_in_from_brain_bit_out(buf_aes, &buf_brain[brain_out_offset]);
     /* init_buffer_aes_p2p_1x1(buf_aes, gold_aes, 0); */

    // input_bytes = val_num * sizeof(token_t);

    ((struct aes_cxx_catapult_access *)cfg_p2p_4x1[4].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_4x1[4].esp_desc)->batch = 4 * N;

    print_brain_bit_config(&cfg_p2p_4x1[0]);
    print_brain_bit_config(&cfg_p2p_4x1[1]);
    print_brain_bit_config(&cfg_p2p_4x1[2]);
    print_brain_bit_config(&cfg_p2p_4x1[3]);
    print_aes_config(&cfg_p2p_4x1[4]);

    printf("\n  ** START **\n");
    esp_run(cfg_p2p_4x1, 5);
    total_time = max(max(cfg_p2p_4x1[0].hw_ns,
                         max(cfg_p2p_4x1[1].hw_ns,
                             max(cfg_p2p_4x1[2].hw_ns,
                                 cfg_p2p_4x1[3].hw_ns))),
                     cfg_p2p_4x1[4].hw_ns);
    printf("\n  ** DONE **\n");

    // errors = validate_buffer_aes(buf_aes, &buf_aes[aes_out_offset], gold_aes, 0);
    printf("INFO: total errors %u\n", errors);

    if (!errors)
        printf("  + TEST PASS\n");
    else
        printf("  + TEST FAIL\n");

    printf("==================================================\n");
    printf("==   Finish ** brain_bit + aes (p2p_4x1) **     ==\n");
    printf("==     N: %d\tTotal time: %llu\n", N, total_time);
    printf("==================================================\n");

    fprintf(log_file, "-- p2p_4x1[0] time: %llu\n", cfg_p2p_4x1[0].hw_ns);
    fprintf(log_file, "-- p2p_4x1[1] time: %llu\n", cfg_p2p_4x1[1].hw_ns);
    fprintf(log_file, "-- p2p_4x1[2] time: %llu\n", cfg_p2p_4x1[2].hw_ns);
    fprintf(log_file, "-- p2p_4x1[3] time: %llu\n", cfg_p2p_4x1[3].hw_ns);
    fprintf(log_file, "-- p2p_4x1[4] time: %llu\n", cfg_p2p_4x1[4].hw_ns);
    fprintf(log_file, "-- p2p_4x1 N: %d time:\t%llu\n", N, total_time);

    esp_free(buf_brain);
    esp_free(buf_aes);

    return total_time;
}

int run_both_p2p_4x4(int N)
{
    printf("==================================================\n");
    printf("==   Start ** brain_bit + aes (p2p_4x4) **      ==\n");
    printf("==================================================\n");

    unsigned long long total_time = 0;

    unsigned errors = 0;

    token_t *buf_brain;
    /* token_t *gold_brain; */

    token_t *buf_aes;
    /* token_t *gold_aes; */

    // set brain_bit config registers
    // key_length = 128;
    // key_batch = 3;
    // key_num = 1;
    // val_num = 0; // key_num * 4;
    // tot_iter = N;
    key_batch = N*1.1;
    key_num = N;
    val_num = 0;
    tot_iter = 1;


    init_parameters_brain_bit();
    // init_parameters_aes_from_brain(val_num);

    buf_brain = (token_t *)esp_alloc(brain_size);
    buf_aes = (token_t *)esp_alloc(brain_size);
    cfg_p2p_4x4[0].hw_buf = buf_brain;
    cfg_p2p_4x4[1].hw_buf = buf_brain;
    cfg_p2p_4x4[2].hw_buf = buf_brain;
    cfg_p2p_4x4[3].hw_buf = buf_brain;
    cfg_p2p_4x4[4].hw_buf = buf_aes;
    cfg_p2p_4x4[5].hw_buf = buf_aes;
    cfg_p2p_4x4[6].hw_buf = buf_aes;
    cfg_p2p_4x4[7].hw_buf = buf_aes;
    memset(buf_brain, 0, brain_size);
    memset(buf_aes, 0, brain_size);

    /* gold_brain = malloc(brain_out_size * N); */
    /* gold_aes = malloc(aes_out_size); */

    // init_buffer_brain_bit(buf_brain, gold_brain);

    unsigned avg_u = *avg_ptr;
    unsigned std_u = *std_ptr;
    unsigned R_u = *R_ptr;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[0].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[0].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[0].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[0].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[0].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[0].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[0].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[0].esp_desc)->tot_iter = tot_iter;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[1].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[1].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[1].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[1].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[1].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[1].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[1].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[1].esp_desc)->tot_iter = tot_iter;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[2].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[2].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[2].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[2].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[2].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[2].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[2].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[2].esp_desc)->tot_iter = tot_iter;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[3].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[3].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[3].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[3].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[3].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[3].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[3].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_4x4[3].esp_desc)->tot_iter = tot_iter;

    set_aes_in_from_brain_bit_out(buf_aes, &buf_brain[brain_out_offset]);
    /* init_buffer_aes_p2p_1x1(buf_aes, gold_aes, 0); */

    // input_bytes = val_num * sizeof(token_t);

    ((struct aes_cxx_catapult_access *)cfg_p2p_4x4[4].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_4x4[4].esp_desc)->batch = N;
    ((struct aes_cxx_catapult_access *)cfg_p2p_4x4[5].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_4x4[5].esp_desc)->batch = N;
    ((struct aes_cxx_catapult_access *)cfg_p2p_4x4[6].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_4x4[6].esp_desc)->batch = N;
    ((struct aes_cxx_catapult_access *)cfg_p2p_4x4[7].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_4x4[7].esp_desc)->batch = N;

    print_brain_bit_config(&cfg_p2p_4x4[0]);
    print_brain_bit_config(&cfg_p2p_4x4[1]);
    print_brain_bit_config(&cfg_p2p_4x4[2]);
    print_brain_bit_config(&cfg_p2p_4x4[3]);
    print_aes_config(&cfg_p2p_4x4[4]);
    print_aes_config(&cfg_p2p_4x4[5]);
    print_aes_config(&cfg_p2p_4x4[6]);
    print_aes_config(&cfg_p2p_4x4[7]);

    printf("\n  ** START **\n");
    esp_run(cfg_p2p_4x4, 8);
    total_time = max(max(cfg_p2p_4x4[0].hw_ns,
                         max(cfg_p2p_4x4[1].hw_ns,
                             max(cfg_p2p_4x4[2].hw_ns,
                                 cfg_p2p_4x4[3].hw_ns))),
                     max(cfg_p2p_4x4[4].hw_ns,
                         max(cfg_p2p_4x4[5].hw_ns,
                             max(cfg_p2p_4x4[6].hw_ns,
                                 cfg_p2p_4x4[7].hw_ns))));
    printf("\n  ** DONE **\n");

    // errors = validate_buffer_aes(buf_aes, &buf_aes[aes_out_offset], gold_aes, 0);
    printf("INFO: total errors %u\n", errors);

    if (!errors)
        printf("  + TEST PASS\n");
    else
        printf("  + TEST FAIL\n");

    printf("==================================================\n");
    printf("==   Finish ** brain_bit + aes (p2p_4x4) **     ==\n");
    printf("==     N: %d\tTotal time: %llu\n", N, total_time);
    printf("==================================================\n");

    fprintf(log_file, "-- p2p_4x4[0] time: %llu\n", cfg_p2p_4x4[0].hw_ns);
    fprintf(log_file, "-- p2p_4x4[1] time: %llu\n", cfg_p2p_4x4[1].hw_ns);
    fprintf(log_file, "-- p2p_4x4[2] time: %llu\n", cfg_p2p_4x4[2].hw_ns);
    fprintf(log_file, "-- p2p_4x4[3] time: %llu\n", cfg_p2p_4x4[3].hw_ns);
    fprintf(log_file, "-- p2p_4x4[4] time: %llu\n", cfg_p2p_4x4[4].hw_ns);
    fprintf(log_file, "-- p2p_4x4[5] time: %llu\n", cfg_p2p_4x4[5].hw_ns);
    fprintf(log_file, "-- p2p_4x4[6] time: %llu\n", cfg_p2p_4x4[6].hw_ns);
    fprintf(log_file, "-- p2p_4x4[7] time: %llu\n", cfg_p2p_4x4[7].hw_ns);
    fprintf(log_file, "-- p2p_4x4 N: %d time:\t%llu\n", N, total_time);

    esp_free(buf_brain);
    esp_free(buf_aes);

    return total_time;
}

int run_both_p2p_2x3(int N)
{
    printf("==================================================\n");
    printf("==   Start ** brain_bit + aes (p2p_2x3) **      ==\n");
    printf("==================================================\n");

    unsigned long long total_time = 0;

    unsigned errors = 0;

    token_t *buf_brain;
    /* token_t *gold_brain; */

    token_t *buf_aes;
    /* token_t *gold_aes; */

    // set brain_bit config registers
    // key_length = 128;
    // key_batch = 3;
    // key_num = 1;
    // val_num = 0; // key_num * 4;
    // tot_iter = N;
    key_batch = N*1.1;
    key_num = N;
    val_num = 0;
    tot_iter = 1;


    init_parameters_brain_bit();
    // init_parameters_aes_from_brain(val_num);

    buf_brain = (token_t *)esp_alloc(brain_size);
    buf_aes = (token_t *)esp_alloc(brain_size);
    cfg_p2p_2x3[0].hw_buf = buf_brain;
    cfg_p2p_2x3[1].hw_buf = buf_brain;
    cfg_p2p_2x3[2].hw_buf = buf_aes;
    cfg_p2p_2x3[3].hw_buf = buf_aes;
    cfg_p2p_2x3[4].hw_buf = buf_aes;
    memset(buf_brain, 0, brain_size);
    memset(buf_aes, 0, brain_size);

    /* gold_brain = malloc(brain_out_size * tot_iter); */
    /* gold_aes = malloc(aes_out_size * tot_iter); */

    // init_buffer_brain_bit(buf_brain, gold_brain);
    unsigned avg_u = *avg_ptr;
    unsigned std_u = *std_ptr;
    unsigned R_u = *R_ptr;

    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[0].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[0].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[0].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[0].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[0].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[0].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[0].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[0].esp_desc)->tot_iter = tot_iter;

    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[1].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[1].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[1].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[1].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[1].esp_desc)->key_batch = 2 * key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[1].esp_desc)->key_num = 2 * key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[1].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_2x3[1].esp_desc)->tot_iter = tot_iter;

    set_aes_in_from_brain_bit_out(buf_aes, &buf_brain[brain_out_offset]);
    /* init_buffer_aes_p2p_1x1(buf_aes, gold_aes, 0); */

    // input_bytes = val_num * sizeof(token_t);

    ((struct aes_cxx_catapult_access *)cfg_p2p_2x3[2].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_2x3[2].esp_desc)->batch = N;

    ((struct aes_cxx_catapult_access *)cfg_p2p_2x3[3].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_2x3[3].esp_desc)->batch = N;

    ((struct aes_cxx_catapult_access *)cfg_p2p_2x3[4].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_2x3[4].esp_desc)->batch = N;

    print_brain_bit_config(&cfg_p2p_2x3[0]);
    print_brain_bit_config(&cfg_p2p_2x3[1]);
    print_aes_config(&cfg_p2p_2x3[2]);
    print_aes_config(&cfg_p2p_2x3[3]);
    print_aes_config(&cfg_p2p_2x3[4]);

    printf("\n  ** START **\n");
    esp_run(cfg_p2p_2x3, 5);
    total_time = max(max(max(max(cfg_p2p_2x3[0].hw_ns, cfg_p2p_2x3[1].hw_ns), cfg_p2p_2x3[2].hw_ns), cfg_p2p_2x3[3].hw_ns), cfg_p2p_2x3[4].hw_ns);
    printf("\n  ** DONE **\n");

    // errors = validate_buffer_aes(buf_aes, &buf_aes[aes_out_offset], gold_aes, 0);
    printf("INFO: total errors %u\n", errors);

    if (!errors)
        printf("  + TEST PASS\n");
    else
        printf("  + TEST FAIL\n");

    printf("==================================================\n");
    printf("==   Finish ** brain_bit + aes (p2p_2x3) **     ==\n");
    printf("==     N: %d\tTotal time: %llu\n", N, total_time);
    printf("==================================================\n");

    fprintf(log_file, "-- p2p_2x3[0] time: %llu\n", cfg_p2p_2x3[0].hw_ns);
    fprintf(log_file, "-- p2p_2x3[1] time: %llu\n", cfg_p2p_2x3[1].hw_ns);
    fprintf(log_file, "-- p2p_2x3[2] time: %llu\n", cfg_p2p_2x3[2].hw_ns);
    fprintf(log_file, "-- p2p_2x3[3] time: %llu\n", cfg_p2p_2x3[3].hw_ns);
    fprintf(log_file, "-- p2p_2x3[4] time: %llu\n", cfg_p2p_2x3[4].hw_ns);
    fprintf(log_file, "-- p2p_2x3 N: %d time:\t%llu\n", N, total_time);

    esp_free(buf_brain);
    esp_free(buf_aes);

    return total_time;
}

int run_both_p2p_3x2(int N)
{
    printf("==================================================\n");
    printf("==   Start ** brain_bit + aes (p2p_3x2) **      ==\n");
    printf("==================================================\n");

    unsigned long long total_time = 0;

    unsigned errors = 0;

    token_t *buf_brain;
    /* token_t *gold_brain; */

    token_t *buf_aes;
    /* token_t *gold_aes; */

    // set brain_bit config registers
    key_length = 128;
    key_batch = 3;
    key_num = 1;
    val_num = key_num * 4;
    tot_iter = N;

    init_parameters_brain_bit();
    // init_parameters_aes_from_brain(val_num);

    buf_brain = (token_t *)esp_alloc(brain_size);
    buf_aes = (token_t *)esp_alloc(brain_size);
    cfg_p2p_3x2[0].hw_buf = buf_brain;
    cfg_p2p_3x2[1].hw_buf = buf_brain;
    cfg_p2p_3x2[2].hw_buf = buf_brain;
    cfg_p2p_3x2[3].hw_buf = buf_aes;
    cfg_p2p_3x2[4].hw_buf = buf_aes;
    memset(buf_brain, 0, brain_size);
    memset(buf_aes, 0, brain_size);

    /* gold_brain = malloc(brain_out_size * N); */
    /* gold_aes = malloc(aes_out_size); */

    // init_buffer_brain_bit(buf_brain, gold_brain);

    unsigned avg_u = *avg_ptr;
    unsigned std_u = *std_ptr;
    unsigned R_u = *R_ptr;

    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->tot_iter = tot_iter;

    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->tot_iter = tot_iter;

    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->avg = avg_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->std = std_u;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->R = R_u;

    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->key_length = key_length;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->key_batch = key_batch;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->key_num = key_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->val_num = val_num;
    ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->tot_iter = tot_iter;

    set_aes_in_from_brain_bit_out(buf_aes, &buf_brain[brain_out_offset]);
    /* init_buffer_aes_p2p_1x1(buf_aes, gold_aes, 0); */

    input_bytes = val_num * sizeof(token_t);

    ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[3].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[3].esp_desc)->batch = 2 * N;
    ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[4].esp_desc)->input_bytes = input_bytes;
    ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[4].esp_desc)->batch = N;

    printf("\n====== %s ====== config registers: \n", cfg_p2p_3x2[0].devname);
    printf("  .avg          = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->avg);
    printf("  .key_length   = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->key_length);
    printf("  .std          = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->std);
    printf("  .R            = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->R);
    printf("  .L            = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->L);
    printf("  .key_batch    = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->key_batch);
    printf("  .key_num      = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->key_num);
    printf("  .val_num      = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->val_num);
    printf("  .tot_iter 	= %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[0].esp_desc)->tot_iter);

    printf("\n====== %s ====== config registers: \n", cfg_p2p_3x2[1].devname);
    printf("  .avg          = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->avg);
    printf("  .key_length   = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->key_length);
    printf("  .std          = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->std);
    printf("  .R            = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->R);
    printf("  .L            = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->L);
    printf("  .key_batch    = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->key_batch);
    printf("  .key_num      = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->key_num);
    printf("  .val_num      = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->val_num);
    printf("  .tot_iter 	= %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[1].esp_desc)->tot_iter);

    printf("\n====== %s ====== config registers: \n", cfg_p2p_3x2[2].devname);
    printf("  .avg          = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->avg);
    printf("  .key_length   = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->key_length);
    printf("  .std          = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->std);
    printf("  .R            = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->R);
    printf("  .L            = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->L);
    printf("  .key_batch    = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->key_batch);
    printf("  .key_num      = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->key_num);
    printf("  .val_num      = %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->val_num);
    printf("  .tot_iter 	= %d\n", ((struct brain_bit_vivado_access *)cfg_p2p_3x2[2].esp_desc)->tot_iter);

    printf("\n====== %s ====== config registers: \n", cfg_p2p_3x2[3].devname);
    printf("  .oper_mode    = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[3].esp_desc)->oper_mode);
    printf("  .encryption   = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[3].esp_desc)->encryption);
    printf("  .key_bytes    = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[3].esp_desc)->key_bytes);
    printf("  .input_bytes  = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[3].esp_desc)->input_bytes);
    printf("  .iv_bytes     = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[3].esp_desc)->iv_bytes);
    printf("  .aad_bytes    = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[3].esp_desc)->aad_bytes);
    printf("  .tag_bytes    = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[3].esp_desc)->tag_bytes);
    printf("  .batch        = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[3].esp_desc)->batch);

    printf("\n====== %s ====== config registers: \n", cfg_p2p_3x2[4].devname);
    printf("  .oper_mode    = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[4].esp_desc)->oper_mode);
    printf("  .encryption   = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[4].esp_desc)->encryption);
    printf("  .key_bytes    = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[4].esp_desc)->key_bytes);
    printf("  .input_bytes  = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[4].esp_desc)->input_bytes);
    printf("  .iv_bytes     = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[4].esp_desc)->iv_bytes);
    printf("  .aad_bytes    = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[4].esp_desc)->aad_bytes);
    printf("  .tag_bytes    = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[4].esp_desc)->tag_bytes);
    printf("  .batch        = %d\n", ((struct aes_cxx_catapult_access *)cfg_p2p_3x2[4].esp_desc)->batch);

    printf("\n  ** START **\n");
    esp_run(cfg_p2p_3x2, 5);
    total_time = max(cfg_p2p_3x2[0].hw_ns,
                     max(cfg_p2p_3x2[1].hw_ns,
                         max(cfg_p2p_3x2[2].hw_ns,
                             max(cfg_p2p_3x2[3].hw_ns,
                                 cfg_p2p_3x2[4].hw_ns))));
    printf("\n  ** DONE **\n");

    /* errors = validate_buffer_aes(buf_aes, &buf_aes[aes_out_offset], gold_aes, 0); */
    printf("INFO: total errors %u\n", errors);

    if (!errors)
        printf("  + TEST PASS\n");
    else
        printf("  + TEST FAIL\n");

    printf("==================================================\n");
    printf("==   Finish ** brain_bit + aes (p2p_3x2) **     ==\n");
    printf("==     N: %d\tTotal time: %llu\n", N, total_time);
    printf("==================================================\n");

    esp_free(buf_brain);

    return errors;
}

int wrap_1x1_mem_avg10(void)
{
    int i, j;
    int N[6] = {1, 10, 50, 100, 500, 1000};
    unsigned long long avg_time;

    for (j = 0; j < 6; j++)
    {
        avg_time = 0;
        for (i = 0; i < 10; i++)
        {
            printf("%d\n", i);
            avg_time += run_both_mem_1x1(N[j]);
        }
        fprintf(log_0309, "-- mem_1x1 key_length = %d, input_bytes = %d, N: %d time:\t%llu\n",
                key_length, input_bytes, N[j], avg_time / 10);
    }

    return 0;
}

int wrap_1x1_avg10(void)
{
    int i, j;
    int N[6] = {1, 10, 50, 100, 500, 1000};
    unsigned long long avg_time;

    for (j = 0; j < 6; j++)
    {
        avg_time = 0;
        for (i = 0; i < 10; i++)
        {
            printf("%d\n", i);
            avg_time += run_both_p2p_1x1(N[j]);
        }
        fprintf(log_0309, "-- p2p_1x1 key_length = %d, input_bytes = %d, N: %d time:\t%llu\n",
                key_length, input_bytes, N[j], avg_time / 10);
    }

    return 0;
}

int wrap_2x1_avg10(void)
{
    int i, j;
    int N[6] = {1, 10, 50, 100, 500, 1000};
    unsigned long long avg_time;

    for (j = 0; j < 6; j++)
    {
        avg_time = 0;
        for (i = 0; i < 10; i++)
        {
            printf("%d\n", i);
            avg_time += run_both_p2p_2x1(N[j]);
        }
        fprintf(log_0309, "-- p2p_2x1 key_length = %d, input_bytes = %d, N: %d time:\t%llu\n",
                key_length, input_bytes, N[j], avg_time / 10);
    }

    return 0;
}

int wrap_3x1_avg10(void)
{
    int i, j;
    int N[6] = {1, 10, 50, 100, 500, 1000};
    unsigned long long avg_time;

    for (j = 0; j < 6; j++)
    {
        avg_time = 0;
        for (i = 0; i < 10; i++)
        {
            printf("%d\n", i);
            avg_time += run_both_p2p_3x1(N[j]);
        }
        fprintf(log_0309, "-- p2p_3x1 key_length = %d, input_bytes = %d, N: %d time:\t%llu\n",
                key_length, input_bytes, N[j], avg_time / 10);
    }

    return 0;
}

int wrap_4x1_avg10(void)
{
    int i, j;
    int N[6] = {1, 10, 50, 100, 500, 1000};
    unsigned long long avg_time;

    for (j = 0; j < 6; j++)
    {
        avg_time = 0;
        for (i = 0; i < 10; i++)
        {
            printf("%d\n", i);
            avg_time += run_both_p2p_4x1(N[j]);
        }
        fprintf(log_0309, "-- p2p_4x1 key_length = %d, input_bytes = %d, N: %d time:\t%llu\n",
                key_length, input_bytes, N[j], avg_time / 10);
    }

    return 0;
}

int wrap_4x4_avg10(void)
{
    int i, j;
    int N[6] = {1, 10, 50, 100, 500, 1000};
    unsigned long long avg_time;

    for (j = 0; j < 6; j++)
    {
        avg_time = 0;
        for (i = 0; i < 10; i++)
        {
            printf("%d\n", i);
            avg_time += run_both_p2p_4x4(N[j]);
        }
        fprintf(log_0309, "-- p2p_4x4 key_length = %d, input_bytes = %d, N: %d time:\t%llu\n",
                key_length, input_bytes, N[j], avg_time / 10);
    }

    return 0;
}

int wrap_1x2_avg10(void)
{
    int i, j;
    int N[6] = {1, 10, 50, 100, 500, 1000};
    unsigned long long avg_time;

    for (j = 0; j < 6; j++)
    {
        avg_time = 0;
        for (i = 0; i < 10; i++)
        {
            printf("%d\n", i);
            avg_time += run_both_p2p_1x2(N[j]);
        }
        fprintf(log_0309, "-- p2p_1x2 key_length = %d, input_bytes = %d, N: %d time:\t%llu\n",
                key_length, input_bytes, N[j], avg_time / 10);
    }

    return 0;
}

int wrap_2x3_avg10(void)
{
    int i, j;
    int N[6] = {1, 10, 50, 100, 500, 1000};
    unsigned long long avg_time;

    for (j = 0; j < 6; j++)
    {
        avg_time = 0;
        for (i = 0; i < 10; i++)
        {
            printf("%d\n", i);
            avg_time += run_both_p2p_2x3(N[j]);
        }
        fprintf(log_0309, "-- p2p_2x3 key_length = %d, input_bytes = %d, N: %d time:\t%llu\n",
                key_length, input_bytes, N[j], avg_time / 10);
    }

    return 0;
}

int main(int argc, char **argv)
{
    log_file = fopen("log.txt", "w");
    log_0309 = fopen("log_0309.txt", "w");

    int total_errors = -7;
    /* int errors_0 = 0; */
    /* int errors_1 = 0; */
    /* int errors_2 = 0; */
    /* int errors_3 = 0; */
    /* int errors_4 = 0; */

    /* int avg_time; */
    //int i, j;
    //int N[6] = {1, 10, 50, 100, 500, 1000};

    // errors_0 = run_brain_bit_only();

    // errors_1 = run_aes_only(N_BATCH);

    // errors_2 = run_both_mem_1x1(1);
    // errors_2 = run_both_mem_1x1(10);
    // errors_2 = run_both_mem_1x1(100);

    // key_length = 256;
    // input_bytes = 16;
    // errors_3 = run_both_p2p_1x1(1);
    // errors_3 = run_both_p2p_2x1(1);
    // errors_3 = run_both_p2p_3x1(1);
    // errors_3 = run_both_p2p_4x1(1);

    // -- 1 brain 1 aes mem
    key_length = 256;
    input_bytes = 16;
    wrap_1x1_mem_avg10();

    key_length = 512;
    input_bytes = 48;
    wrap_1x1_mem_avg10();

    key_length = 1024;
    input_bytes = 112;
    wrap_1x1_mem_avg10();

    // key_length = 2048;
    // input_bytes = 240;
    // wrap_1x1_mem_avg10();

    // // -- 1 brain 1 aes
    // key_length = 256;
    // input_bytes = 16;
    // wrap_1x1_avg10();

    // key_length = 512;
    // input_bytes = 48;
    // wrap_1x1_avg10();

    // key_length = 1024;
    // input_bytes = 112;
    // wrap_1x1_avg10();

    // // -- 2 brain 1 aes
    // key_length = 256;
    // input_bytes = 16;
    // wrap_2x1_avg10();

    // key_length = 512;
    // input_bytes = 48;
    // wrap_2x1_avg10();

    // key_length = 1024;
    // input_bytes = 112;
    // wrap_2x1_avg10();

    // // -- 3 brain 1 aes
    // key_length = 256;
    // input_bytes = 16;
    // wrap_3x1_avg10();

    // key_length = 512;
    // input_bytes = 48;
    // wrap_3x1_avg10();

    // key_length = 1024;
    // input_bytes = 112;
    // wrap_3x1_avg10();

    // // -- 4 brain 1 aes
    // key_length = 256;
    // input_bytes = 16;
    // wrap_4x1_avg10();

    // key_length = 512;
    // input_bytes = 48;
    // wrap_4x1_avg10();

    // key_length = 1024;
    // input_bytes = 112;
    // wrap_4x1_avg10();

    // // -- 4 brain 4 aes
    // key_length = 256;
    // input_bytes = 16;
    // wrap_4x4_avg10();

    // key_length = 512;
    // input_bytes = 48;
    // wrap_4x4_avg10();

    // key_length = 1024;
    // input_bytes = 112;
    // wrap_4x4_avg10();

    // // -- 1 brain 2 aes
    // key_length = 256;
    // input_bytes = 16;
    // wrap_1x2_avg10();

    // key_length = 512;
    // input_bytes = 48;
    // wrap_1x2_avg10();

    // key_length = 1024;
    // input_bytes = 112;
    // wrap_1x2_avg10();

    // -- 2 brain 3 aes
    // key_length = 256;
    // input_bytes = 16;
    // wrap_2x3_avg10();

    // key_length = 512;
    // input_bytes = 48;
    // wrap_2x3_avg10();

    // key_length = 1024;
    // input_bytes = 112;
    // wrap_2x3_avg10();

    // errors_4 = run_both_p2p_1x2();

    // total_errors = errors_0 + errors_1 + errors_2 + errors_3 + errors_4;
    // printf("\n---> summary: errors_0 = %d, errors_1 = %d, errors_2 = %d, errors_3 = %d\n",
    //        errors_0, errors_1, errors_2, errors_3);

    fclose(log_file);
    fclose(log_0309);

    printf("0317.04\n");

    return total_errors;
}

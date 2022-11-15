// Copyright (c) 2011-2021 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "libesp.h"
#include "cfg.h"
#include "data.h"

#include <sys/time.h>

#define N_BATCH 7 
static unsigned in_words_adj;
static unsigned out_words_adj;

static unsigned in_len;
static unsigned out_len;
//static unsigned in_size;
//static unsigned out_size;
static unsigned out_offset;
static unsigned size_bytes;


//static unsigned tag_bytes;
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
static unsigned in_size;
static unsigned out_size;
/* static unsigned aad_size; */
/* static unsigned tag_size; */

//input_size for each iteration
const unsigned array_size[12] = {16, 32, 48, 64, 80, 96, 112, 128, 144};

/* User-defined code */
static int validate_buffer(token_t *out, token_t *gold, unsigned indx)
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

        printf("INFO: [%u] @%p %x (%x) %s\n", j, out + j, out_data, gold_data, ((out_data != gold_data)?" !!!":""));

        if (out_data != gold_data)
        {
            errors++;
        //    printf("INFO: [%u] @%p %x (%x) %s\n", j, out + j, out_data, gold_data, ((out_data != gold_data)?" !!!":""));
        }
    }

    printf("INFO: total errors %u\n", errors);

    return errors;
}


/* User-defined code */
static void init_buffer(token_t *in, token_t * gold, token_t *out, unsigned indx)
{
    int i;
    int j;

    printf("INFO: raw_encrypt_plaintext_words[%u] %u\n", indx, ecb_raw_encrypt_plaintext_words[indx]);
    printf("INFO: ecb_raw_encrypt_ciphertext_words[%u] %u\n", indx, ecb_raw_encrypt_ciphertext_words[indx]);
    
    for (j = 0; j < ecb_raw_encrypt_key_words[indx]; j++) {
        in[j] = ecb_raw_encrypt_key[indx][j];
        printf("INFO: raw_encrypt_key[%u][%u] | %x\n", indx, j, in[j]);
    }
    
    for (i = 0; i < ecb_raw_encrypt_plaintext_words[indx]; i++, j++)
    {
        in[j] = ecb_raw_encrypt_plaintext[indx][i];
        printf("INFO: raw_encrypt_plaintext[%u][%u] | inputs[%u]@%p %x\n", indx, i, j, in + j, in[j]);
    }
    
    for (i = 0; i < ecb_raw_encrypt_plaintext_words[indx]; i++, j++)
    {
        //in[j] = ecb_raw_encrypt_plaintext[indx][i];
        printf("INFO: raw_encrypty_output[%u][%u] | outputs[%u]@%p %x\n", indx, i, i, out + i, *(out + i));
    }

    for (j = 0; j < ecb_raw_encrypt_ciphertext_words[indx]; j++) {
        gold[j] = ecb_raw_encrypt_ciphertext[indx][j];
        printf("INFO: raw_encrypt_ciphertext[%u][%u] %x\n", indx, j, gold[j]);
    }

}


/* User-defined code */
static void init_parameters(unsigned indx)
{
    in_words = ecb_raw_encrypt_plaintext_words[indx] + 
                       ecb_raw_encrypt_key_words[indx]; 
    out_words = ecb_raw_encrypt_plaintext_words[indx] ;

    if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0) {
        in_words_adj = in_words;
        out_words_adj = out_words;
    } else {
        in_words_adj = round_up(in_words, DMA_WORD_PER_BEAT(sizeof(token_t)));
        out_words_adj = round_up(out_words, DMA_WORD_PER_BEAT(sizeof(token_t)));
    }
    
    in_len = in_words_adj;
    printf("%s: in_len = %u\n", __func__, in_len);
    out_len =  out_words_adj;
    printf("%s: out_len = %u\n", __func__, out_len);
    in_size = in_len * sizeof(token_t);
    printf("%s: in_size = %u\n", __func__, in_size);
    out_size = out_len * sizeof(token_t);
    printf("%s: out_size = %u\n", __func__, out_size);
    out_offset = in_len;
    printf("%s: out_offset = %u\n", __func__, out_offset);
    size_bytes = (out_offset * sizeof(token_t)) + out_size;
    printf("%s: size = %u\n", __func__, size_bytes);
/*
    in_size = ecb_raw_encrypt_plaintext_bytes [indx];
    key_size = ecb_raw_encrypt_key_bytes[indx];
    out_size = in_size;
    out_offset = ecb_raw_encrypt_plaintext_words[indx] + ecb_raw_encrypt_key_words[indx];
    size_bytes = in_size + key_size + out_size;
    printf("%s: in_size = %u out_size =%u  key_size = %u  tot_size =%u\n", __func__, in_size, out_size, key_size, size_bytes);
 */
 }

int main(int argc, char **argv)
{
    unsigned errors_0 = 0;
    unsigned errors_1 = 0;

    token_t *gold;
    token_t *buf_0;
    /* token_t *buf_1; */

    printf("INFO:   sizeof(token_t) = %lu\n", sizeof(token_t));

    for (unsigned i = 1; i < N_BATCH; i++) {
        init_parameters(i);
        /* const uint32_t new_in_size = in_size; */
        /* volatile uint32_t new_in_size_reg = *(volatile uint32_t *) &new_in_size; // = in_size - 16; */

        buf_0 = (token_t *) esp_alloc(size_bytes);
        cfg_000[0].hw_buf = buf_0;

        memset(buf_0, 0, size_bytes);


        gold = malloc(out_size);

        init_buffer(buf_0, gold, &buf_0[out_offset], i);

        printf("\n====== %s ======\n\n", cfg_000[0].devname);
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


        ((struct aes_cxx_catapult_access*) cfg_000[0].esp_desc)->input_bytes = array_size[i]; //new_in_size_reg; //in_size - 16;

        printf("oper_mode %u\n", ((struct aes_cxx_catapult_access*) cfg_000[0].esp_desc)->oper_mode);
        printf("encryption %u\n", ((struct aes_cxx_catapult_access*) cfg_000[0].esp_desc)->encryption);
        printf("key_bytes  %u\n", ((struct aes_cxx_catapult_access*) cfg_000[0].esp_desc)->key_bytes);
        printf("input_bytes %u\n", ((struct aes_cxx_catapult_access*) cfg_000[0].esp_desc)->input_bytes);
        printf("iv_bytes %u\n", ((struct aes_cxx_catapult_access*) cfg_000[0].esp_desc)->iv_bytes);
        printf("aad_bytes %u\n", ((struct aes_cxx_catapult_access*) cfg_000[0].esp_desc)->aad_bytes);
        printf("tag_bytes %u\n", ((struct aes_cxx_catapult_access*) cfg_000[0].esp_desc)->tag_bytes);

        struct timeval  hw_begin_0, hw_end_0;
        gettimeofday(&hw_begin_0, NULL);
        esp_run(cfg_000, NACC);
        gettimeofday(&hw_end_0, NULL);

        printf("\n  ** DONE **\n");

        errors_0 = validate_buffer(&buf_0[out_offset], gold, i);

        free(gold);
        esp_free(buf_0);

        if (!errors_0)
            printf("  + TEST PASS\n");
        else
            printf("  + TEST FAIL\n");

        printf("\n====== %s ======\n\n", cfg_000[0].devname);
    }

    return (errors_0 + errors_1);
}


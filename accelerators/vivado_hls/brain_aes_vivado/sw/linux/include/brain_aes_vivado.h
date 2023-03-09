// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef _BRAIN_AES_VIVADO_H_
#define _BRAIN_AES_VIVADO_H_

#ifdef __KERNEL__
#include <linux/ioctl.h>
#include <linux/types.h>
#else
#include <sys/ioctl.h>
#include <stdint.h>
#ifndef __user
#define __user
#endif
#endif /* __KERNEL__ */

#include <esp.h>
#include <esp_accelerator.h>

#include <stdio.h>

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

FILE *log_file;

struct brain_bit_vivado_access
{
	struct esp_access esp;
	/* <<--regs-->> */
	unsigned avg;
	unsigned key_length;
	unsigned std;
	unsigned R;
	unsigned L;
	unsigned key_batch;
	unsigned key_num;
	unsigned val_num;
	unsigned tot_iter;
	unsigned src_offset;
	unsigned dst_offset;
};

struct aes_cxx_catapult_access
{
	struct esp_access esp;
	/* <<--regs-->> */
	unsigned oper_mode;
	unsigned encryption;
	unsigned key_bytes;
	unsigned input_bytes;
	unsigned iv_bytes;
	unsigned aad_bytes;
	unsigned tag_bytes;
	unsigned batch;
	unsigned src_offset;
	unsigned dst_offset;
};

#define BRAIN_BIT_VIVADO_IOC_ACCESS _IOW('S', 0, struct brain_bit_vivado_access)
#define AES_CXX_CATAPULT_IOC_ACCESS _IOW('S', 0, struct aes_cxx_catapult_access)

typedef int32_t token_t;
#define DATA_BITWIDTH 32
//#define N_TESTS 4

/* <<--params-def-->> */
// brain_bit:
#define AVG 3.0677295382679177
#define KEY_LENGTH 128
#define STD 38.626628825256695
#define R_val 1.5
#define L_val 1500
#define KEY_BATCH 20
#define KEY_NUM 15
#define VAL_NUM 1
#define TOT_ITER 1
// aes:
#define OPER_MODE 1
#define ENCRYPTION 1
#define KEY_BYTES 16
#define INPUT_BYTES 8
#define IV_BYTES 0
#define AAD_BYTES 0
#define TAG_BYTES 0
#define BATCH 1

/* <<--params-->> */
// brain_bit:
const float avg = AVG;
unsigned *avg_ptr = (unsigned *)&avg;
int32_t key_length = KEY_LENGTH;
const float std = STD;
unsigned *std_ptr = (unsigned *)&std;
const float R = R_val;
unsigned *R_ptr = (unsigned *)&R;
const int32_t L = L_val;
int32_t key_batch = KEY_BATCH;
int32_t key_num = KEY_NUM;
int32_t val_num = VAL_NUM;
int32_t tot_iter = TOT_ITER;
const float Rs = R * std;
// aes:
unsigned oper_mode = OPER_MODE;
unsigned encryption = ENCRYPTION;
unsigned key_bytes = KEY_BYTES;
unsigned input_bytes = INPUT_BYTES;
unsigned iv_bytes = IV_BYTES;
unsigned aad_bytes = AAD_BYTES;
unsigned tag_bytes = TAG_BYTES;
unsigned batch = BATCH;

#define NACC 1


// brain_bit:
 unsigned brain_in_words_adj;
 unsigned brain_out_words_adj;
 unsigned brain_in_len;
 unsigned brain_out_len;
 unsigned brain_in_size;
 unsigned brain_out_size;
 unsigned brain_out_offset;
 unsigned brain_size;
 int brain_key_counter;

// aes:
#define N_BATCH 9
 unsigned aes_in_words_adj;
 unsigned aes_out_words_adj;
 unsigned aes_in_len;
 unsigned aes_out_len;
 unsigned aes_in_size;
 unsigned aes_out_size;
 unsigned aes_out_offset;
 unsigned aes_size_bytes;

//  unsigned tag_bytes;
/*  unsigned aad_bytes; */
/*  unsigned in_bytes; */
/*  unsigned out_bytes; */
/*  unsigned iv_bytes; */
/*  unsigned key_bytes; */
/*  unsigned encryption; */
/*  unsigned oper_mode; */

/*  unsigned key_words; */
/*  unsigned iv_words; */
 unsigned aes_in_words;
 unsigned aes_out_words;
/*  unsigned aad_words; */
/*  unsigned tag_words; */

/*  unsigned key_size; */
/*  unsigned iv_size; */
//  unsigned in_size;
//  unsigned out_size;
/*  unsigned aad_size; */
/*  unsigned tag_size; */


// for run brain_bit only
void init_parameters_brain_bit(void);
void init_buffer_brain_bit(token_t *in, token_t *gold);
int validate_buffer_brain_bit(token_t *out, token_t *gold);

// for run aes only
void init_parameters_aes(unsigned indx);
void init_buffer_aes(token_t *in, token_t *gold, token_t *out, unsigned indx);
int validate_buffer_aes(token_t *in, token_t *out, token_t *gold, unsigned indx);


// for run both mem 1x1
void init_parameters_aes_from_brain(int val_n);
void set_aes_in_from_brain_bit_out(token_t *in_aes, token_t *out_brain);
// void init_buffer_aes_from_brain(token_t *in, token_t *aes_key, token_t *aes_val, token_t *out, unsigned indx);


void init_buffer_aes_p2p_1_1(token_t *in, token_t *gold, unsigned indx);


#endif /* _BRAIN_AES_VIVADO_H_ */

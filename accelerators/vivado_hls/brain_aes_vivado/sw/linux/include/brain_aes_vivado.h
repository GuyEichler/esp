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

struct brain_bit_vivado_access {
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

struct aes_cxx_catapult_access {
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


#define BRAIN_BIT_VIVADO_IOC_ACCESS	_IOW ('S', 0, struct brain_bit_vivado_access)
#define AES_CXX_CATAPULT_IOC_ACCESS	_IOW ('S', 0, struct aes_cxx_catapult_access)

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

#endif /* _BRAIN_AES_VIVADO_H_ */

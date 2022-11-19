// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef __CFG_P2P_H__
#define __CFG_P2P_H__

#include "libesp.h"
#include "brain_aes_vivado.h"

typedef int32_t token_t;

// /* <<--params-def-->> */
// // brain_bit:
// #define AVG 3.0677295382679177
// #define KEY_LENGTH 128
// #define STD 38.626628825256695
// #define R_val 1.5
// #define L_val 1500
// #define KEY_BATCH 20
// #define KEY_NUM 15
// #define VAL_NUM 1
// #define TOT_ITER 1
// // aes:
// #define OPER_MODE 1
// #define ENCRYPTION 1
// #define KEY_BYTES 16
// #define INPUT_BYTES 8
// #define IV_BYTES 0
// #define AAD_BYTES 0
// #define TAG_BYTES 0
// #define BATCH 1

// /* <<--params-->> */
// // brain_bit:
// const float avg = AVG;
// unsigned *avg_ptr = (unsigned *)&avg;
// const int32_t key_length = KEY_LENGTH;
// const float std = STD;
// unsigned *std_ptr = (unsigned *)&std;
// const float R = R_val;
// unsigned *R_ptr = (unsigned *)&R;
// const int32_t L = L_val;
// const int32_t key_batch = KEY_BATCH;
// const int32_t key_num = KEY_NUM;
// const int32_t val_num = VAL_NUM;
// const int32_t tot_iter = TOT_ITER;
// const float Rs = R * std;
// // aes:
// const unsigned oper_mode = OPER_MODE;
// const unsigned encryption = ENCRYPTION;
// const unsigned key_bytes = KEY_BYTES;
// const unsigned input_bytes = INPUT_BYTES;
// const unsigned iv_bytes = IV_BYTES;
// const unsigned aad_bytes = AAD_BYTES;
// const unsigned tag_bytes = TAG_BYTES;
// const unsigned batch = BATCH;

// #define NACC 1

struct brain_bit_vivado_access brain_bit_cfg_p2p_000[] = {
	{
		/* <<--descriptor-->> */
		.avg = AVG,
		.key_length = KEY_LENGTH,
		.std = STD,
		.R = R,
		.L = L,
		.key_batch = KEY_BATCH,
		.key_num = KEY_NUM,
		.val_num = VAL_NUM,
		.tot_iter = TOT_ITER,
		.src_offset = 0,
		.dst_offset = 0,
		.esp.coherence = ACC_COH_NONE,
		.esp.p2p_store = 1,
		.esp.p2p_nsrcs = 0,
		.esp.p2p_srcs = {"", "", "", ""},
	}};

struct aes_cxx_catapult_access aes_cxx_cfg_p2p_000[] = {
	{
		/* <<--descriptor-->> */
		.oper_mode = OPER_MODE,
		.encryption = ENCRYPTION,
		.key_bytes = KEY_BYTES,
		.input_bytes = INPUT_BYTES,
		.iv_bytes = IV_BYTES,
		.aad_bytes = AAD_BYTES,
		.tag_bytes = TAG_BYTES,
		.batch = BATCH,
		.src_offset = 0,
		.dst_offset = 0,
		.esp.coherence = ACC_COH_NONE,
		.esp.p2p_store = 0,
		.esp.p2p_nsrcs = 1,
		.esp.p2p_srcs = {"brain_bit_vivado.0", "", "", ""},
	}};

esp_thread_info_t cfg_p2p_000[] = {
	{
		.run = true,
		.devname = "brain_bit_vivado.0",
		.ioctl_req = BRAIN_BIT_VIVADO_IOC_ACCESS,
		.esp_desc = &(brain_bit_cfg_000[0].esp),
	},
	{
		.run = true,
		.devname = "aes_cxx_catapult.0",
		.ioctl_req = AES_CXX_CATAPULT_IOC_ACCESS,
		.esp_desc = &(aes_cxx_cfg_000[0].esp),
	}};

#endif /* __CFG_P2P_H__ */

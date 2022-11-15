// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef __ESP_CFG_000_H__
#define __ESP_CFG_000_H__

#include "libesp.h"
#include "brain_aes_vivado.h"

typedef int32_t token_t;

/* <<--params-def-->> */
#define AVG 3.0677295382679177
#define KEY_LENGTH 128
#define STD 38.626628825256695
#define R_val 1.5
#define L_val 1500
#define KEY_BATCH 20
#define KEY_NUM 15
#define VAL_NUM 0

/* <<--params-->> */
const float avg = AVG;
unsigned* avg_ptr = (unsigned*)&avg;
const int32_t key_length = KEY_LENGTH;
const float std = STD;
unsigned* std_ptr = (unsigned*)&std;
const float R = R_val;
unsigned* R_ptr = (unsigned*)&R;
const int32_t L = L_val;
const int32_t key_batch = KEY_BATCH;
const int32_t key_num = KEY_NUM;
const int32_t val_num = VAL_NUM;
const float Rs = R * std;

#define NACC 1

struct brain_bit_vivado_access brain_bit_cfg_000[] = {
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
		.src_offset = 0,
		.dst_offset = 0,
		.esp.coherence = ACC_COH_NONE,
		.esp.p2p_store = 0,
		.esp.p2p_nsrcs = 0,
		.esp.p2p_srcs = {"", "", "", ""},
	}
};

esp_thread_info_t cfg_000[] = {
	{
		.run = true,
		.devname = "brain_bit_vivado.0",
		.ioctl_req = BRAIN_BIT_VIVADO_IOC_ACCESS,
		.esp_desc = &(brain_bit_cfg_000[0].esp),
	}
};

#endif /* __ESP_CFG_000_H__ */

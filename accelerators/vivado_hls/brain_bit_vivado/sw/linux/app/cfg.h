// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef __ESP_CFG_000_H__
#define __ESP_CFG_000_H__

#include "libesp.h"
#include "brain_bit_vivado.h"

typedef int32_t token_t;

/* <<--params-def-->> */
#define AVG 1
#define KEY_LENGTH 128
#define STD 1
#define R 1
#define L 1
#define KEY_BATCH 1

/* <<--params-->> */
const int32_t avg = AVG;
const int32_t key_length = KEY_LENGTH;
const int32_t std = STD;
const int32_t R = R;
const int32_t L = L;
const int32_t key_batch = KEY_BATCH;

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

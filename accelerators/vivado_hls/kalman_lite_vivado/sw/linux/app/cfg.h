// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef __ESP_CFG_000_H__
#define __ESP_CFG_000_H__

#include "libesp.h"
#include "kalman_lite_vivado.h"

typedef float token_t;

#define STATES 6
#define NEURONS 164
#define TIME_STAMPS 100
#define CHUNKS 1
#define BATCHES TIME_STAMPS / CHUNKS

/* <<--params-def-->> */
#define ITER BATCHES
#define X_DIM STATES
#define Z_DIM NEURONS
#define INV_NUM 2
#define INV_RESET 3

/* <<--params-->> */
const int32_t inv_reset = INV_RESET;
const int32_t inv_num = INV_NUM;
const int32_t chunks = CHUNKS;
const int32_t iter = BATCHES;
const int32_t x_dim = X_DIM;
const int32_t z_dim = Z_DIM;

#define NACC 1

struct kalman_lite_vivado_access kalman_lite_cfg_000[] = {
	{
		/* <<--descriptor-->> */
		.inv_reset = INV_RESET,
		.inv_num = INV_NUM,
		.chunks = CHUNKS,
		.iter = BATCHES,
		.x_dim = X_DIM,
		.z_dim = Z_DIM,
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
		.devname = "kalman_lite_vivado.0",
		.ioctl_req = KALMAN_LITE_VIVADO_IOC_ACCESS,
		.esp_desc = &(kalman_lite_cfg_000[0].esp),
	}
};

#endif /* __ESP_CFG_000_H__ */

// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef __ESP_CFG_000_H__
#define __ESP_CFG_000_H__

#include "libesp.h"
#include "kalman_gauss_vivado.h"

typedef float token_t;

#define STATES 2
#define NEURONS 1
#define TIME_STAMPS 2
#define CHUNKS 1
#define BATCHES TIME_STAMPS / CHUNKS

/* <<--params-def-->> */
#define ITER BATCHES
#define X_DIM STATES
#define Z_DIM NEURONS
#define INV_NUM 0
#define INV_RESET 0

/* <<--params-->> */
int32_t inv_reset = INV_RESET;
int32_t inv_num = INV_NUM;
int32_t chunks = CHUNKS;
int32_t iter = BATCHES;
int32_t x_dim = X_DIM;
int32_t z_dim = Z_DIM;

#define NACC 1

struct kalman_gauss_vivado_access kalman_gauss_cfg_000[] = {
	{
		/* <<--descriptor-->> */
		.inv_reset = INV_RESET,
		.inv_num = INV_NUM,
		.chunks = CHUNKS,
		.iter = ITER,
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
		.devname = "kalman_gauss_vivado.0",
		.ioctl_req = KALMAN_GAUSS_VIVADO_IOC_ACCESS,
		.esp_desc = &(kalman_gauss_cfg_000[0].esp),
	}
};

#endif /* __ESP_CFG_000_H__ */

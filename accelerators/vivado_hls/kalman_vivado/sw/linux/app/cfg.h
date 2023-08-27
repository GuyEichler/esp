// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef __ESP_CFG_000_H__
#define __ESP_CFG_000_H__

#include "libesp.h"
#include "kalman_vivado.h"

typedef int32_t token_t;

/* <<--params-def-->> */
#define ITER 1
#define X_DIM 1
#define Z_DIM 1

/* <<--params-->> */
const int32_t iter = ITER;
const int32_t x_dim = X_DIM;
const int32_t z_dim = Z_DIM;

#define NACC 1

struct kalman_vivado_access kalman_cfg_000[] = {
	{
		/* <<--descriptor-->> */
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
		.devname = "kalman_vivado.0",
		.ioctl_req = KALMAN_VIVADO_IOC_ACCESS,
		.esp_desc = &(kalman_cfg_000[0].esp),
	}
};

#endif /* __ESP_CFG_000_H__ */

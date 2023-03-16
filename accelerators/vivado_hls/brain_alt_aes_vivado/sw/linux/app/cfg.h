// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef __ESP_CFG_000_H__
#define __ESP_CFG_000_H__

#include "libesp.h"
#include "brain_alt_aes_vivado.h"

struct brain_bit_alt_vivado_access brain_bit_alt_cfg_000[] = {
	{
		/* <<--descriptor-->> */
		.avg = AVG,
		.key_length = KEY_LENGTH,
		.std = STD,
		.R = R,
		.d = D,
		.h = H,
		.key_batch = KEY_BATCH,
		.key_num = KEY_NUM,
		.val_num = VAL_NUM,
		.tot_iter = TOT_ITER,
		.src_offset = 0,
		.dst_offset = 0,
		.esp.coherence = ACC_COH_NONE,
		.esp.p2p_store = 0,
		.esp.p2p_nsrcs = 0,
		.esp.p2p_srcs = {"", "", "", ""},
	}};

struct aes_cxx_catapult_access aes_cxx_cfg_000[] = {
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
		.esp.p2p_nsrcs = 0,
		.esp.p2p_srcs = {"", "", "", ""},
	}};

esp_thread_info_t cfg_brain_bit_alt_000[] = {
	{
		.run = true,
		.devname = "brain_bit_alt_vivado.0",
		.ioctl_req = BRAIN_BIT_ALT_VIVADO_IOC_ACCESS,
		.esp_desc = &(brain_bit_alt_cfg_000[0].esp),
	}};

esp_thread_info_t cfg_aes_000[] = {
	{
		.run = true,
		.devname = "aes_cxx_catapult.0",
		.ioctl_req = AES_CXX_CATAPULT_IOC_ACCESS,
		.esp_desc = &(aes_cxx_cfg_000[0].esp),
	}};

#endif /* __ESP_CFG_000_H__ */

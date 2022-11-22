// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef __CFG_P2P_1x4_H__
#define __CFG_P2P_1x4_H__

#include "libesp.h"
#include "brain_aes_vivado.h"

struct brain_bit_vivado_access brain_bit_cfg_p2p_1x4[] = {
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

struct aes_cxx_catapult_access aes_cxx_cfg_p2p_1x4[] = {
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
	},
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
	},
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
	},
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

esp_thread_info_t cfg_p2p_1x4[] = {
	{
		.run = true,
		.devname = "brain_bit_vivado.0",
		.ioctl_req = BRAIN_BIT_VIVADO_IOC_ACCESS,
		.esp_desc = &(brain_bit_cfg_p2p_1x4[0].esp),
	},
	{
		.run = true,
		.devname = "aes_cxx_catapult.0",
		.ioctl_req = AES_CXX_CATAPULT_IOC_ACCESS,
		.esp_desc = &(aes_cxx_cfg_p2p_1x4[0].esp),
	},
	{
		.run = true,
		.devname = "aes_cxx_catapult.1",
		.ioctl_req = AES_CXX_CATAPULT_IOC_ACCESS,
		.esp_desc = &(aes_cxx_cfg_p2p_1x4[1].esp),
	},
	{
		.run = true,
		.devname = "aes_cxx_catapult.2",
		.ioctl_req = AES_CXX_CATAPULT_IOC_ACCESS,
		.esp_desc = &(aes_cxx_cfg_p2p_1x4[2].esp),
	},
	{
		.run = true,
		.devname = "aes_cxx_catapult.3",
		.ioctl_req = AES_CXX_CATAPULT_IOC_ACCESS,
		.esp_desc = &(aes_cxx_cfg_p2p_1x4[3].esp),
	}
	};

#endif /* __CFG_P2P_1x4_H__ */

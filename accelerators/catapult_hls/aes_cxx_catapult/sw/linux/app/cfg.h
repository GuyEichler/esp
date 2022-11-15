// Copyright (c) 2011-2021 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef __ESP_CFG_000_H__
#define __ESP_CFG_000_H__

#include "libesp.h"
#include "aes_cxx_catapult.h"

typedef int32_t token_t;

/* <<--params-def-->> */
#define OPER_MODE 1
#define ENCRYPTION 1
#define KEY_BYTES 16
#define INPUT_BYTES  8
#define IV_BYTES  0
#define AAD_BYTES 0
#define TAG_BYTES 0
#define BATCH 1



/* <<--params-->> */
const unsigned oper_mode    = OPER_MODE;
const unsigned encryption   = ENCRYPTION;
const unsigned key_bytes    = KEY_BYTES;
const unsigned input_bytes  = INPUT_BYTES;
const unsigned iv_bytes     = IV_BYTES;
const unsigned aad_bytes    = AAD_BYTES;
const unsigned tag_bytes    = TAG_BYTES;
const unsigned batch    = BATCH;

/* << --hardcoded -->> */

#define NACC 1


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
     }
};


esp_thread_info_t cfg_000[] = {
    {
        .run = true,
        .devname = "aes_cxx_catapult.0",
        .ioctl_req = AES_CXX_CATAPULT_IOC_ACCESS,
        .esp_desc = &(aes_cxx_cfg_000[0].esp),
    }
};

#endif /* __ESP_CFG_000_H__ */

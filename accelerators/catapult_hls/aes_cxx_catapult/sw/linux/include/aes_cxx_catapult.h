// Copyright (c) 2011-2021 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef _AES_CXX_CATAPULT_H_
#define _AES_CXX_CATAPULT_H_

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

#define AES_CXX_CATAPULT_IOC_ACCESS	_IOW ('S', 0, struct aes_cxx_catapult_access)

#endif /* _AES_CXX_CATAPULT_H_ */

// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef _BRAIN_BIT_ALT_VIVADO_H_
#define _BRAIN_BIT_ALT_VIVADO_H_

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

struct brain_bit_alt_vivado_access {
	struct esp_access esp;
	/* <<--regs-->> */
	unsigned avg;
	unsigned key_length;
	unsigned std;
	unsigned R;
	/* unsigned L; */
	unsigned key_batch;
	unsigned key_num;
	unsigned val_num;
	unsigned tot_iter;
	unsigned d;
	unsigned h;
	unsigned src_offset;
	unsigned dst_offset;
};

#define BRAIN_BIT_ALT_VIVADO_IOC_ACCESS	_IOW ('S', 0, struct brain_bit_alt_vivado_access)

#endif /* _BRAIN_BIT_ALT_VIVADO_H_ */

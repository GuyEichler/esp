// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef _KALMAN_FX_VIVADO_H_
#define _KALMAN_FX_VIVADO_H_

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

struct kalman_fx_vivado_access {
	struct esp_access esp;
	/* <<--regs-->> */
	unsigned inv_reset;
	unsigned inv_num;
	unsigned chunks;
	unsigned iter;
	unsigned x_dim;
	unsigned z_dim;
	unsigned src_offset;
	unsigned dst_offset;
};

#define KALMAN_FX_VIVADO_IOC_ACCESS	_IOW ('S', 0, struct kalman_fx_vivado_access)

#endif /* _KALMAN_FX_VIVADO_H_ */

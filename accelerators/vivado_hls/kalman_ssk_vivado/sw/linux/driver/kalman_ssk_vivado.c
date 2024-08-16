// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include <linux/of_device.h>
#include <linux/mm.h>

#include <asm/io.h>

#include <esp_accelerator.h>
#include <esp.h>

#include "kalman_ssk_vivado.h"

#define DRV_NAME	"kalman_ssk_vivado"

/* <<--regs-->> */
#define KALMAN_SSK_INV_RESET_REG 0x54
#define KALMAN_SSK_INV_NUM_REG 0x50
#define KALMAN_SSK_CHUNKS_REG 0x4C
#define KALMAN_SSK_ITER_REG 0x48
#define KALMAN_SSK_X_DIM_REG 0x44
#define KALMAN_SSK_Z_DIM_REG 0x40

struct kalman_ssk_vivado_device {
	struct esp_device esp;
};

static struct esp_driver kalman_ssk_driver;

static struct of_device_id kalman_ssk_device_ids[] = {
	{
		.name = "SLD_KALMAN_SSK_VIVADO",
	},
	{
		.name = "eb_249",
	},
	{
		.compatible = "sld,kalman_ssk_vivado",
	},
	{ },
};

static int kalman_ssk_devs;

static inline struct kalman_ssk_vivado_device *to_kalman_ssk(struct esp_device *esp)
{
	return container_of(esp, struct kalman_ssk_vivado_device, esp);
}

static void kalman_ssk_prep_xfer(struct esp_device *esp, void *arg)
{
	struct kalman_ssk_vivado_access *a = arg;

	/* <<--regs-config-->> */
	iowrite32be(a->inv_reset, esp->iomem + KALMAN_SSK_INV_RESET_REG);
	iowrite32be(a->inv_num, esp->iomem + KALMAN_SSK_INV_NUM_REG);
	iowrite32be(a->chunks, esp->iomem + KALMAN_SSK_CHUNKS_REG);
	iowrite32be(a->iter, esp->iomem + KALMAN_SSK_ITER_REG);
	iowrite32be(a->x_dim, esp->iomem + KALMAN_SSK_X_DIM_REG);
	iowrite32be(a->z_dim, esp->iomem + KALMAN_SSK_Z_DIM_REG);
	iowrite32be(a->src_offset, esp->iomem + SRC_OFFSET_REG);
	iowrite32be(a->dst_offset, esp->iomem + DST_OFFSET_REG);

}

static bool kalman_ssk_xfer_input_ok(struct esp_device *esp, void *arg)
{
	/* struct kalman_ssk_vivado_device *kalman_ssk = to_kalman_ssk(esp); */
	/* struct kalman_ssk_vivado_access *a = arg; */

	return true;
}

static int kalman_ssk_probe(struct platform_device *pdev)
{
	struct kalman_ssk_vivado_device *kalman_ssk;
	struct esp_device *esp;
	int rc;

	kalman_ssk = kzalloc(sizeof(*kalman_ssk), GFP_KERNEL);
	if (kalman_ssk == NULL)
		return -ENOMEM;
	esp = &kalman_ssk->esp;
	esp->module = THIS_MODULE;
	esp->number = kalman_ssk_devs;
	esp->driver = &kalman_ssk_driver;
	rc = esp_device_register(esp, pdev);
	if (rc)
		goto err;

	kalman_ssk_devs++;
	return 0;
 err:
	kfree(kalman_ssk);
	return rc;
}

static int __exit kalman_ssk_remove(struct platform_device *pdev)
{
	struct esp_device *esp = platform_get_drvdata(pdev);
	struct kalman_ssk_vivado_device *kalman_ssk = to_kalman_ssk(esp);

	esp_device_unregister(esp);
	kfree(kalman_ssk);
	return 0;
}

static struct esp_driver kalman_ssk_driver = {
	.plat = {
		.probe		= kalman_ssk_probe,
		.remove		= kalman_ssk_remove,
		.driver		= {
			.name = DRV_NAME,
			.owner = THIS_MODULE,
			.of_match_table = kalman_ssk_device_ids,
		},
	},
	.xfer_input_ok	= kalman_ssk_xfer_input_ok,
	.prep_xfer	= kalman_ssk_prep_xfer,
	.ioctl_cm	= KALMAN_SSK_VIVADO_IOC_ACCESS,
	.arg_size	= sizeof(struct kalman_ssk_vivado_access),
};

static int __init kalman_ssk_init(void)
{
	return esp_driver_register(&kalman_ssk_driver);
}

static void __exit kalman_ssk_exit(void)
{
	esp_driver_unregister(&kalman_ssk_driver);
}

module_init(kalman_ssk_init)
module_exit(kalman_ssk_exit)

MODULE_DEVICE_TABLE(of, kalman_ssk_device_ids);

MODULE_AUTHOR("Emilio G. Cota <cota@braap.org>");
MODULE_LICENSE("GPL");
MODULE_DESCRIPTION("kalman_ssk_vivado driver");

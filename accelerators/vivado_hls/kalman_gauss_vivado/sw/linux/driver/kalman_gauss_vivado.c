// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include <linux/of_device.h>
#include <linux/mm.h>

#include <asm/io.h>

#include <esp_accelerator.h>
#include <esp.h>

#include "kalman_gauss_vivado.h"

#define DRV_NAME	"kalman_gauss_vivado"

/* <<--regs-->> */
#define KALMAN_GAUSS_INV_RESET_REG 0x54
#define KALMAN_GAUSS_INV_NUM_REG 0x50
#define KALMAN_GAUSS_CHUNKS_REG 0x4C
#define KALMAN_GAUSS_ITER_REG 0x48
#define KALMAN_GAUSS_X_DIM_REG 0x44
#define KALMAN_GAUSS_Z_DIM_REG 0x40

struct kalman_gauss_vivado_device {
	struct esp_device esp;
};

static struct esp_driver kalman_gauss_driver;

static struct of_device_id kalman_gauss_device_ids[] = {
	{
		.name = "SLD_KALMAN_GAUSS_VIVADO",
	},
	{
		.name = "eb_051",
	},
	{
		.compatible = "sld,kalman_gauss_vivado",
	},
	{ },
};

static int kalman_gauss_devs;

static inline struct kalman_gauss_vivado_device *to_kalman_gauss(struct esp_device *esp)
{
	return container_of(esp, struct kalman_gauss_vivado_device, esp);
}

static void kalman_gauss_prep_xfer(struct esp_device *esp, void *arg)
{
	struct kalman_gauss_vivado_access *a = arg;

	/* <<--regs-config-->> */
	iowrite32be(a->inv_reset, esp->iomem + KALMAN_GAUSS_INV_RESET_REG);
	iowrite32be(a->inv_num, esp->iomem + KALMAN_GAUSS_INV_NUM_REG);
	iowrite32be(a->chunks, esp->iomem + KALMAN_GAUSS_CHUNKS_REG);
	iowrite32be(a->iter, esp->iomem + KALMAN_GAUSS_ITER_REG);
	iowrite32be(a->x_dim, esp->iomem + KALMAN_GAUSS_X_DIM_REG);
	iowrite32be(a->z_dim, esp->iomem + KALMAN_GAUSS_Z_DIM_REG);
	iowrite32be(a->src_offset, esp->iomem + SRC_OFFSET_REG);
	iowrite32be(a->dst_offset, esp->iomem + DST_OFFSET_REG);

}

static bool kalman_gauss_xfer_input_ok(struct esp_device *esp, void *arg)
{
	/* struct kalman_gauss_vivado_device *kalman_gauss = to_kalman_gauss(esp); */
	/* struct kalman_gauss_vivado_access *a = arg; */

	return true;
}

static int kalman_gauss_probe(struct platform_device *pdev)
{
	struct kalman_gauss_vivado_device *kalman_gauss;
	struct esp_device *esp;
	int rc;

	kalman_gauss = kzalloc(sizeof(*kalman_gauss), GFP_KERNEL);
	if (kalman_gauss == NULL)
		return -ENOMEM;
	esp = &kalman_gauss->esp;
	esp->module = THIS_MODULE;
	esp->number = kalman_gauss_devs;
	esp->driver = &kalman_gauss_driver;
	rc = esp_device_register(esp, pdev);
	if (rc)
		goto err;

	kalman_gauss_devs++;
	return 0;
 err:
	kfree(kalman_gauss);
	return rc;
}

static int __exit kalman_gauss_remove(struct platform_device *pdev)
{
	struct esp_device *esp = platform_get_drvdata(pdev);
	struct kalman_gauss_vivado_device *kalman_gauss = to_kalman_gauss(esp);

	esp_device_unregister(esp);
	kfree(kalman_gauss);
	return 0;
}

static struct esp_driver kalman_gauss_driver = {
	.plat = {
		.probe		= kalman_gauss_probe,
		.remove		= kalman_gauss_remove,
		.driver		= {
			.name = DRV_NAME,
			.owner = THIS_MODULE,
			.of_match_table = kalman_gauss_device_ids,
		},
	},
	.xfer_input_ok	= kalman_gauss_xfer_input_ok,
	.prep_xfer	= kalman_gauss_prep_xfer,
	.ioctl_cm	= KALMAN_GAUSS_VIVADO_IOC_ACCESS,
	.arg_size	= sizeof(struct kalman_gauss_vivado_access),
};

static int __init kalman_gauss_init(void)
{
	return esp_driver_register(&kalman_gauss_driver);
}

static void __exit kalman_gauss_exit(void)
{
	esp_driver_unregister(&kalman_gauss_driver);
}

module_init(kalman_gauss_init)
module_exit(kalman_gauss_exit)

MODULE_DEVICE_TABLE(of, kalman_gauss_device_ids);

MODULE_AUTHOR("Emilio G. Cota <cota@braap.org>");
MODULE_LICENSE("GPL");
MODULE_DESCRIPTION("kalman_gauss_vivado driver");

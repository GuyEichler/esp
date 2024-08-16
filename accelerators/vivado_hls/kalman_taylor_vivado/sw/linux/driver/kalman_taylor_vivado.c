// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include <linux/of_device.h>
#include <linux/mm.h>

#include <asm/io.h>

#include <esp_accelerator.h>
#include <esp.h>

#include "kalman_taylor_vivado.h"

#define DRV_NAME	"kalman_taylor_vivado"

/* <<--regs-->> */
#define KALMAN_TAYLOR_INV_RESET_REG 0x54
#define KALMAN_TAYLOR_INV_NUM_REG 0x50
#define KALMAN_TAYLOR_CHUNKS_REG 0x4C
#define KALMAN_TAYLOR_ITER_REG 0x48
#define KALMAN_TAYLOR_X_DIM_REG 0x44
#define KALMAN_TAYLOR_Z_DIM_REG 0x40

struct kalman_taylor_vivado_device {
	struct esp_device esp;
};

static struct esp_driver kalman_taylor_driver;

static struct of_device_id kalman_taylor_device_ids[] = {
	{
		.name = "SLD_KALMAN_TAYLOR_VIVADO",
	},
	{
		.name = "eb_350",
	},
	{
		.compatible = "sld,kalman_taylor_vivado",
	},
	{ },
};

static int kalman_taylor_devs;

static inline struct kalman_taylor_vivado_device *to_kalman_taylor(struct esp_device *esp)
{
	return container_of(esp, struct kalman_taylor_vivado_device, esp);
}

static void kalman_taylor_prep_xfer(struct esp_device *esp, void *arg)
{
	struct kalman_taylor_vivado_access *a = arg;

	/* <<--regs-config-->> */
	iowrite32be(a->inv_reset, esp->iomem + KALMAN_TAYLOR_INV_RESET_REG);
	iowrite32be(a->inv_num, esp->iomem + KALMAN_TAYLOR_INV_NUM_REG);
	iowrite32be(a->chunks, esp->iomem + KALMAN_TAYLOR_CHUNKS_REG);
	iowrite32be(a->iter, esp->iomem + KALMAN_TAYLOR_ITER_REG);
	iowrite32be(a->x_dim, esp->iomem + KALMAN_TAYLOR_X_DIM_REG);
	iowrite32be(a->z_dim, esp->iomem + KALMAN_TAYLOR_Z_DIM_REG);
	iowrite32be(a->src_offset, esp->iomem + SRC_OFFSET_REG);
	iowrite32be(a->dst_offset, esp->iomem + DST_OFFSET_REG);

}

static bool kalman_taylor_xfer_input_ok(struct esp_device *esp, void *arg)
{
	/* struct kalman_taylor_vivado_device *kalman_taylor = to_kalman_taylor(esp); */
	/* struct kalman_taylor_vivado_access *a = arg; */

	return true;
}

static int kalman_taylor_probe(struct platform_device *pdev)
{
	struct kalman_taylor_vivado_device *kalman_taylor;
	struct esp_device *esp;
	int rc;

	kalman_taylor = kzalloc(sizeof(*kalman_taylor), GFP_KERNEL);
	if (kalman_taylor == NULL)
		return -ENOMEM;
	esp = &kalman_taylor->esp;
	esp->module = THIS_MODULE;
	esp->number = kalman_taylor_devs;
	esp->driver = &kalman_taylor_driver;
	rc = esp_device_register(esp, pdev);
	if (rc)
		goto err;

	kalman_taylor_devs++;
	return 0;
 err:
	kfree(kalman_taylor);
	return rc;
}

static int __exit kalman_taylor_remove(struct platform_device *pdev)
{
	struct esp_device *esp = platform_get_drvdata(pdev);
	struct kalman_taylor_vivado_device *kalman_taylor = to_kalman_taylor(esp);

	esp_device_unregister(esp);
	kfree(kalman_taylor);
	return 0;
}

static struct esp_driver kalman_taylor_driver = {
	.plat = {
		.probe		= kalman_taylor_probe,
		.remove		= kalman_taylor_remove,
		.driver		= {
			.name = DRV_NAME,
			.owner = THIS_MODULE,
			.of_match_table = kalman_taylor_device_ids,
		},
	},
	.xfer_input_ok	= kalman_taylor_xfer_input_ok,
	.prep_xfer	= kalman_taylor_prep_xfer,
	.ioctl_cm	= KALMAN_TAYLOR_VIVADO_IOC_ACCESS,
	.arg_size	= sizeof(struct kalman_taylor_vivado_access),
};

static int __init kalman_taylor_init(void)
{
	return esp_driver_register(&kalman_taylor_driver);
}

static void __exit kalman_taylor_exit(void)
{
	esp_driver_unregister(&kalman_taylor_driver);
}

module_init(kalman_taylor_init)
module_exit(kalman_taylor_exit)

MODULE_DEVICE_TABLE(of, kalman_taylor_device_ids);

MODULE_AUTHOR("Emilio G. Cota <cota@braap.org>");
MODULE_LICENSE("GPL");
MODULE_DESCRIPTION("kalman_taylor_vivado driver");

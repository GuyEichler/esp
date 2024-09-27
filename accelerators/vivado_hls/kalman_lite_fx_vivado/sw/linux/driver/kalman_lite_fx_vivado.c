// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include <linux/of_device.h>
#include <linux/mm.h>

#include <asm/io.h>

#include <esp_accelerator.h>
#include <esp.h>

#include "kalman_lite_fx_vivado.h"

#define DRV_NAME	"kalman_lite_fx_vivado"

/* <<--regs-->> */
#define KALMAN_LITE_FX_INV_RESET_REG 0x54
#define KALMAN_LITE_FX_INV_NUM_REG 0x50
#define KALMAN_LITE_FX_CHUNKS_REG 0x4C
#define KALMAN_LITE_FX_ITER_REG 0x48
#define KALMAN_LITE_FX_X_DIM_REG 0x44
#define KALMAN_LITE_FX_Z_DIM_REG 0x40

struct kalman_lite_fx_vivado_device {
	struct esp_device esp;
};

static struct esp_driver kalman_lite_fx_driver;

static struct of_device_id kalman_lite_fx_device_ids[] = {
	{
		.name = "SLD_KALMAN_LITE_FX_VIVADO",
	},
	{
		.name = "eb_252",
	},
	{
		.compatible = "sld,kalman_lite_fx_vivado",
	},
	{ },
};

static int kalman_lite_fx_devs;

static inline struct kalman_lite_fx_vivado_device *to_kalman_lite_fx(struct esp_device *esp)
{
	return container_of(esp, struct kalman_lite_fx_vivado_device, esp);
}

static void kalman_lite_fx_prep_xfer(struct esp_device *esp, void *arg)
{
	struct kalman_lite_fx_vivado_access *a = arg;

	/* <<--regs-config-->> */
	iowrite32be(a->inv_reset, esp->iomem + KALMAN_LITE_FX_INV_RESET_REG);
	iowrite32be(a->inv_num, esp->iomem + KALMAN_LITE_FX_INV_NUM_REG);
	iowrite32be(a->chunks, esp->iomem + KALMAN_LITE_FX_CHUNKS_REG);
	iowrite32be(a->iter, esp->iomem + KALMAN_LITE_FX_ITER_REG);
	iowrite32be(a->x_dim, esp->iomem + KALMAN_LITE_FX_X_DIM_REG);
	iowrite32be(a->z_dim, esp->iomem + KALMAN_LITE_FX_Z_DIM_REG);
	iowrite32be(a->src_offset, esp->iomem + SRC_OFFSET_REG);
	iowrite32be(a->dst_offset, esp->iomem + DST_OFFSET_REG);

}

static bool kalman_lite_fx_xfer_input_ok(struct esp_device *esp, void *arg)
{
	/* struct kalman_lite_fx_vivado_device *kalman_lite_fx = to_kalman_lite_fx(esp); */
	/* struct kalman_lite_fx_vivado_access *a = arg; */

	return true;
}

static int kalman_lite_fx_probe(struct platform_device *pdev)
{
	struct kalman_lite_fx_vivado_device *kalman_lite_fx;
	struct esp_device *esp;
	int rc;

	kalman_lite_fx = kzalloc(sizeof(*kalman_lite_fx), GFP_KERNEL);
	if (kalman_lite_fx == NULL)
		return -ENOMEM;
	esp = &kalman_lite_fx->esp;
	esp->module = THIS_MODULE;
	esp->number = kalman_lite_fx_devs;
	esp->driver = &kalman_lite_fx_driver;
	rc = esp_device_register(esp, pdev);
	if (rc)
		goto err;

	kalman_lite_fx_devs++;
	return 0;
 err:
	kfree(kalman_lite_fx);
	return rc;
}

static int __exit kalman_lite_fx_remove(struct platform_device *pdev)
{
	struct esp_device *esp = platform_get_drvdata(pdev);
	struct kalman_lite_fx_vivado_device *kalman_lite_fx = to_kalman_lite_fx(esp);

	esp_device_unregister(esp);
	kfree(kalman_lite_fx);
	return 0;
}

static struct esp_driver kalman_lite_fx_driver = {
	.plat = {
		.probe		= kalman_lite_fx_probe,
		.remove		= kalman_lite_fx_remove,
		.driver		= {
			.name = DRV_NAME,
			.owner = THIS_MODULE,
			.of_match_table = kalman_lite_fx_device_ids,
		},
	},
	.xfer_input_ok	= kalman_lite_fx_xfer_input_ok,
	.prep_xfer	= kalman_lite_fx_prep_xfer,
	.ioctl_cm	= KALMAN_LITE_FX_VIVADO_IOC_ACCESS,
	.arg_size	= sizeof(struct kalman_lite_fx_vivado_access),
};

static int __init kalman_lite_fx_init(void)
{
	return esp_driver_register(&kalman_lite_fx_driver);
}

static void __exit kalman_lite_fx_exit(void)
{
	esp_driver_unregister(&kalman_lite_fx_driver);
}

module_init(kalman_lite_fx_init)
module_exit(kalman_lite_fx_exit)

MODULE_DEVICE_TABLE(of, kalman_lite_fx_device_ids);

MODULE_AUTHOR("Emilio G. Cota <cota@braap.org>");
MODULE_LICENSE("GPL");
MODULE_DESCRIPTION("kalman_lite_fx_vivado driver");

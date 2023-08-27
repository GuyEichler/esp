// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include <linux/of_device.h>
#include <linux/mm.h>

#include <asm/io.h>

#include <esp_accelerator.h>
#include <esp.h>

#include "kalman_vivado.h"

#define DRV_NAME	"kalman_vivado"

/* <<--regs-->> */
#define KALMAN_ITER_REG 0x48
#define KALMAN_X_DIM_REG 0x44
#define KALMAN_Z_DIM_REG 0x40

struct kalman_vivado_device {
	struct esp_device esp;
};

static struct esp_driver kalman_driver;

static struct of_device_id kalman_device_ids[] = {
	{
		.name = "SLD_KALMAN_VIVADO",
	},
	{
		.name = "eb_549",
	},
	{
		.compatible = "sld,kalman_vivado",
	},
	{ },
};

static int kalman_devs;

static inline struct kalman_vivado_device *to_kalman(struct esp_device *esp)
{
	return container_of(esp, struct kalman_vivado_device, esp);
}

static void kalman_prep_xfer(struct esp_device *esp, void *arg)
{
	struct kalman_vivado_access *a = arg;

	/* <<--regs-config-->> */
	iowrite32be(a->iter, esp->iomem + KALMAN_ITER_REG);
	iowrite32be(a->x_dim, esp->iomem + KALMAN_X_DIM_REG);
	iowrite32be(a->z_dim, esp->iomem + KALMAN_Z_DIM_REG);
	iowrite32be(a->src_offset, esp->iomem + SRC_OFFSET_REG);
	iowrite32be(a->dst_offset, esp->iomem + DST_OFFSET_REG);

}

static bool kalman_xfer_input_ok(struct esp_device *esp, void *arg)
{
	/* struct kalman_vivado_device *kalman = to_kalman(esp); */
	/* struct kalman_vivado_access *a = arg; */

	return true;
}

static int kalman_probe(struct platform_device *pdev)
{
	struct kalman_vivado_device *kalman;
	struct esp_device *esp;
	int rc;

	kalman = kzalloc(sizeof(*kalman), GFP_KERNEL);
	if (kalman == NULL)
		return -ENOMEM;
	esp = &kalman->esp;
	esp->module = THIS_MODULE;
	esp->number = kalman_devs;
	esp->driver = &kalman_driver;
	rc = esp_device_register(esp, pdev);
	if (rc)
		goto err;

	kalman_devs++;
	return 0;
 err:
	kfree(kalman);
	return rc;
}

static int __exit kalman_remove(struct platform_device *pdev)
{
	struct esp_device *esp = platform_get_drvdata(pdev);
	struct kalman_vivado_device *kalman = to_kalman(esp);

	esp_device_unregister(esp);
	kfree(kalman);
	return 0;
}

static struct esp_driver kalman_driver = {
	.plat = {
		.probe		= kalman_probe,
		.remove		= kalman_remove,
		.driver		= {
			.name = DRV_NAME,
			.owner = THIS_MODULE,
			.of_match_table = kalman_device_ids,
		},
	},
	.xfer_input_ok	= kalman_xfer_input_ok,
	.prep_xfer	= kalman_prep_xfer,
	.ioctl_cm	= KALMAN_VIVADO_IOC_ACCESS,
	.arg_size	= sizeof(struct kalman_vivado_access),
};

static int __init kalman_init(void)
{
	return esp_driver_register(&kalman_driver);
}

static void __exit kalman_exit(void)
{
	esp_driver_unregister(&kalman_driver);
}

module_init(kalman_init)
module_exit(kalman_exit)

MODULE_DEVICE_TABLE(of, kalman_device_ids);

MODULE_AUTHOR("Emilio G. Cota <cota@braap.org>");
MODULE_LICENSE("GPL");
MODULE_DESCRIPTION("kalman_vivado driver");

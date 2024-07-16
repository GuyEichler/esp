// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include <linux/of_device.h>
#include <linux/mm.h>

#include <asm/io.h>

#include <esp_accelerator.h>
#include <esp.h>

#include "kalman_qr_vivado.h"

#define DRV_NAME	"kalman_qr_vivado"

/* <<--regs-->> */
#define KALMAN_QR_INV_RESET_REG 0x54
#define KALMAN_QR_INV_NUM_REG 0x50
#define KALMAN_QR_CHUNKS_REG 0x4C
#define KALMAN_QR_ITER_REG 0x48
#define KALMAN_QR_X_DIM_REG 0x44
#define KALMAN_QR_Z_DIM_REG 0x40

struct kalman_qr_vivado_device {
	struct esp_device esp;
};

static struct esp_driver kalman_qr_driver;

static struct of_device_id kalman_qr_device_ids[] = {
	{
		.name = "SLD_KALMAN_QR_VIVADO",
	},
	{
		.name = "eb_069",
	},
	{
		.compatible = "sld,kalman_qr_vivado",
	},
	{ },
};

static int kalman_qr_devs;

static inline struct kalman_qr_vivado_device *to_kalman_qr(struct esp_device *esp)
{
	return container_of(esp, struct kalman_qr_vivado_device, esp);
}

static void kalman_qr_prep_xfer(struct esp_device *esp, void *arg)
{
	struct kalman_qr_vivado_access *a = arg;

	/* <<--regs-config-->> */
	iowrite32be(a->inv_reset, esp->iomem + KALMAN_QR_INV_RESET_REG);
	iowrite32be(a->inv_num, esp->iomem + KALMAN_QR_INV_NUM_REG);
	iowrite32be(a->chunks, esp->iomem + KALMAN_QR_CHUNKS_REG);
	iowrite32be(a->iter, esp->iomem + KALMAN_QR_ITER_REG);
	iowrite32be(a->x_dim, esp->iomem + KALMAN_QR_X_DIM_REG);
	iowrite32be(a->z_dim, esp->iomem + KALMAN_QR_Z_DIM_REG);
	iowrite32be(a->src_offset, esp->iomem + SRC_OFFSET_REG);
	iowrite32be(a->dst_offset, esp->iomem + DST_OFFSET_REG);

}

static bool kalman_qr_xfer_input_ok(struct esp_device *esp, void *arg)
{
	/* struct kalman_qr_vivado_device *kalman_qr = to_kalman_qr(esp); */
	/* struct kalman_qr_vivado_access *a = arg; */

	return true;
}

static int kalman_qr_probe(struct platform_device *pdev)
{
	struct kalman_qr_vivado_device *kalman_qr;
	struct esp_device *esp;
	int rc;

	kalman_qr = kzalloc(sizeof(*kalman_qr), GFP_KERNEL);
	if (kalman_qr == NULL)
		return -ENOMEM;
	esp = &kalman_qr->esp;
	esp->module = THIS_MODULE;
	esp->number = kalman_qr_devs;
	esp->driver = &kalman_qr_driver;
	rc = esp_device_register(esp, pdev);
	if (rc)
		goto err;

	kalman_qr_devs++;
	return 0;
 err:
	kfree(kalman_qr);
	return rc;
}

static int __exit kalman_qr_remove(struct platform_device *pdev)
{
	struct esp_device *esp = platform_get_drvdata(pdev);
	struct kalman_qr_vivado_device *kalman_qr = to_kalman_qr(esp);

	esp_device_unregister(esp);
	kfree(kalman_qr);
	return 0;
}

static struct esp_driver kalman_qr_driver = {
	.plat = {
		.probe		= kalman_qr_probe,
		.remove		= kalman_qr_remove,
		.driver		= {
			.name = DRV_NAME,
			.owner = THIS_MODULE,
			.of_match_table = kalman_qr_device_ids,
		},
	},
	.xfer_input_ok	= kalman_qr_xfer_input_ok,
	.prep_xfer	= kalman_qr_prep_xfer,
	.ioctl_cm	= KALMAN_QR_VIVADO_IOC_ACCESS,
	.arg_size	= sizeof(struct kalman_qr_vivado_access),
};

static int __init kalman_qr_init(void)
{
	return esp_driver_register(&kalman_qr_driver);
}

static void __exit kalman_qr_exit(void)
{
	esp_driver_unregister(&kalman_qr_driver);
}

module_init(kalman_qr_init)
module_exit(kalman_qr_exit)

MODULE_DEVICE_TABLE(of, kalman_qr_device_ids);

MODULE_AUTHOR("Emilio G. Cota <cota@braap.org>");
MODULE_LICENSE("GPL");
MODULE_DESCRIPTION("kalman_qr_vivado driver");

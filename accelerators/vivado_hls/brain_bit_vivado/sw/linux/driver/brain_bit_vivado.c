// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include <linux/of_device.h>
#include <linux/mm.h>

#include <asm/io.h>

#include <esp_accelerator.h>
#include <esp.h>

#include "brain_bit_vivado.h"

#define DRV_NAME	"brain_bit_vivado"

/* <<--regs-->> */
#define BRAIN_BIT_TOT_ITER_REG 0x60
#define BRAIN_BIT_VAL_NUM_REG 0x5c
#define BRAIN_BIT_KEY_NUM_REG 0x58
#define BRAIN_BIT_AVG_REG 0x54
#define BRAIN_BIT_KEY_LENGTH_REG 0x50
#define BRAIN_BIT_STD_REG 0x4c
#define BRAIN_BIT_R_REG 0x48
#define BRAIN_BIT_L_REG 0x44
#define BRAIN_BIT_KEY_BATCH_REG 0x40

struct brain_bit_vivado_device {
	struct esp_device esp;
};

static struct esp_driver brain_bit_driver;

static struct of_device_id brain_bit_device_ids[] = {
	{
		.name = "SLD_BRAIN_BIT_VIVADO",
	},
	{
		.name = "eb_1e4",
	},
	{
		.compatible = "sld,brain_bit_vivado",
	},
	{ },
};

static int brain_bit_devs;

static inline struct brain_bit_vivado_device *to_brain_bit(struct esp_device *esp)
{
	return container_of(esp, struct brain_bit_vivado_device, esp);
}

static void brain_bit_prep_xfer(struct esp_device *esp, void *arg)
{
	struct brain_bit_vivado_access *a = arg;

	/* <<--regs-config-->> */
	iowrite32be(a->avg, esp->iomem + BRAIN_BIT_AVG_REG);
	iowrite32be(a->key_length, esp->iomem + BRAIN_BIT_KEY_LENGTH_REG);
	iowrite32be(a->std, esp->iomem + BRAIN_BIT_STD_REG);
	iowrite32be(a->R, esp->iomem + BRAIN_BIT_R_REG);
	iowrite32be(a->L, esp->iomem + BRAIN_BIT_L_REG);
	iowrite32be(a->key_batch, esp->iomem + BRAIN_BIT_KEY_BATCH_REG);
	iowrite32be(a->key_num, esp->iomem + BRAIN_BIT_KEY_NUM_REG);
	iowrite32be(a->val_num, esp->iomem + BRAIN_BIT_VAL_NUM_REG);
	iowrite32be(a->tot_iter, esp->iomem + BRAIN_BIT_TOT_ITER_REG);
	iowrite32be(a->src_offset, esp->iomem + SRC_OFFSET_REG);
	iowrite32be(a->dst_offset, esp->iomem + DST_OFFSET_REG);

}

static bool brain_bit_xfer_input_ok(struct esp_device *esp, void *arg)
{
	/* struct brain_bit_vivado_device *brain_bit = to_brain_bit(esp); */
	/* struct brain_bit_vivado_access *a = arg; */

	return true;
}

static int brain_bit_probe(struct platform_device *pdev)
{
	struct brain_bit_vivado_device *brain_bit;
	struct esp_device *esp;
	int rc;

	brain_bit = kzalloc(sizeof(*brain_bit), GFP_KERNEL);
	if (brain_bit == NULL)
		return -ENOMEM;
	esp = &brain_bit->esp;
	esp->module = THIS_MODULE;
	esp->number = brain_bit_devs;
	esp->driver = &brain_bit_driver;
	rc = esp_device_register(esp, pdev);
	if (rc)
		goto err;

	brain_bit_devs++;
	return 0;
 err:
	kfree(brain_bit);
	return rc;
}

static int __exit brain_bit_remove(struct platform_device *pdev)
{
	struct esp_device *esp = platform_get_drvdata(pdev);
	struct brain_bit_vivado_device *brain_bit = to_brain_bit(esp);

	esp_device_unregister(esp);
	kfree(brain_bit);
	return 0;
}

static struct esp_driver brain_bit_driver = {
	.plat = {
		.probe		= brain_bit_probe,
		.remove		= brain_bit_remove,
		.driver		= {
			.name = DRV_NAME,
			.owner = THIS_MODULE,
			.of_match_table = brain_bit_device_ids,
		},
	},
	.xfer_input_ok	= brain_bit_xfer_input_ok,
	.prep_xfer	= brain_bit_prep_xfer,
	.ioctl_cm	= BRAIN_BIT_VIVADO_IOC_ACCESS,
	.arg_size	= sizeof(struct brain_bit_vivado_access),
};

static int __init brain_bit_init(void)
{
	return esp_driver_register(&brain_bit_driver);
}

static void __exit brain_bit_exit(void)
{
	esp_driver_unregister(&brain_bit_driver);
}

module_init(brain_bit_init)
module_exit(brain_bit_exit)

MODULE_DEVICE_TABLE(of, brain_bit_device_ids);

MODULE_AUTHOR("Emilio G. Cota <cota@braap.org>");
MODULE_LICENSE("GPL");
MODULE_DESCRIPTION("brain_bit_vivado driver");

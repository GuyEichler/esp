// Copyright (c) 2011-2021 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include <linux/of_device.h>
#include <linux/mm.h>

#include <asm/io.h>

#include <esp_accelerator.h>
#include <esp.h>

#include "aes_cxx_catapult.h"

#define DRV_NAME "aes_cxx_catapult"

/* <<--regs-->> */
#define AES_CXX_OPER_MODE_REG 0x40
#define AES_CXX_ENCRYPTION_REG 0x44
#define AES_CXX_KEY_BYTES_REG 0x48
#define AES_CXX_INPUT_BYTES_REG 0x4C
#define AES_CXX_IV_BYTES_REG 0x50
#define AES_CXX_AAD_BYTES_REG 0x54
#define AES_CXX_TAG_BYTES_REG 0x58

struct aes_cxx_catapult_device {
	struct esp_device esp;
};

static struct esp_driver aes_cxx_driver;

static struct of_device_id aes_cxx_device_ids[] = {
	{
		.name = "SLD_AES_CXX_CATAPULT",
	},
	{
		.name = "eb_089",
	},
	{
		.compatible = "sld,aes_cxx_catapult",
	},
	{ },
};

static int aes_cxx_devs;

static inline struct aes_cxx_catapult_device *to_aes_cxx(struct esp_device *esp)
{
	return container_of(esp, struct aes_cxx_catapult_device, esp);
}

static void aes_cxx_prep_xfer(struct esp_device *esp, void *arg)
{
	struct aes_cxx_catapult_access *a = arg;

	/* <<--regs-config-->> */
	iowrite32be(a->oper_mode, esp->iomem + AES_CXX_OPER_MODE_REG);
	iowrite32be(a->encryption, esp->iomem + AES_CXX_ENCRYPTION_REG);
	iowrite32be(a->key_bytes, esp->iomem + AES_CXX_KEY_BYTES_REG);
	iowrite32be(a->input_bytes, esp->iomem + AES_CXX_INPUT_BYTES_REG);
	iowrite32be(a->iv_bytes, esp->iomem + AES_CXX_IV_BYTES_REG);
	iowrite32be(a->aad_bytes, esp->iomem + AES_CXX_AAD_BYTES_REG);
	iowrite32be(a->tag_bytes, esp->iomem + AES_CXX_TAG_BYTES_REG);
	iowrite32be(a->src_offset, esp->iomem + SRC_OFFSET_REG);
	iowrite32be(a->dst_offset, esp->iomem + DST_OFFSET_REG);

}

static bool aes_cxx_xfer_input_ok(struct esp_device *esp, void *arg)
{
	/* struct softmax_cxx_catapult_device *softmax_cxx = to_softmax_cxx(esp); */
	/* struct softmax_cxx_catapult_access *a = arg; */

	return true;
}

static int aes_cxx_probe(struct platform_device *pdev)
{
	struct aes_cxx_catapult_device *aes_cxx;
	struct esp_device *esp;
	int rc;

	aes_cxx = kzalloc(sizeof(*aes_cxx), GFP_KERNEL);
	if (aes_cxx == NULL)
		return -ENOMEM;
	esp = &aes_cxx->esp;
	esp->module = THIS_MODULE;
	esp->number = aes_cxx_devs;
	esp->driver = &aes_cxx_driver;
	rc = esp_device_register(esp, pdev);
	if (rc)
		goto err;

	aes_cxx_devs++;
	return 0;
 err:
	kfree(aes_cxx);
	return rc;
}

static int __exit aes_cxx_remove(struct platform_device *pdev)
{
	struct esp_device *esp = platform_get_drvdata(pdev);
	struct aes_cxx_catapult_device *aes_cxx = to_aes_cxx(esp);

	esp_device_unregister(esp);
	kfree(aes_cxx);
	return 0;
}

static struct esp_driver aes_cxx_driver = {
	.plat = {
		.probe		= aes_cxx_probe,
		.remove		= aes_cxx_remove,
		.driver		= {
			.name = DRV_NAME,
			.owner = THIS_MODULE,
			.of_match_table = aes_cxx_device_ids,
		},
	},
	.xfer_input_ok	= aes_cxx_xfer_input_ok,
	.prep_xfer	= aes_cxx_prep_xfer,
	.ioctl_cm	= AES_CXX_CATAPULT_IOC_ACCESS,
	.arg_size	= sizeof(struct aes_cxx_catapult_access),
};

static int __init aes_cxx_init(void)
{
	return esp_driver_register(&aes_cxx_driver);
}

static void __exit aes_cxx_exit(void)
{
	esp_driver_unregister(&aes_cxx_driver);
}

module_init(aes_cxx_init)
module_exit(aes_cxx_exit)

MODULE_DEVICE_TABLE(of, aes_cxx_device_ids);

MODULE_AUTHOR("Giuseppe Di Guglielmo <giuseppe@cs.columbia.edu>");
MODULE_LICENSE("GPL");
MODULE_DESCRIPTION("softmax_cxx_catapult driver");

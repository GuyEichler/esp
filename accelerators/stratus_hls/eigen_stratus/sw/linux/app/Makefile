# Copyright (c) 2011-2023 Columbia University, System Level Design Group
# SPDX-License-Identifier: Apache-2.0

EXTRA_OBJS := libcrungemm.so
CFLAGS += -L.
LDFLAGS += -lcrungemm
CXX := riscv64-unknown-linux-gnu-g++

APPNAME := eigen_gemm
include $(DRIVERS)/common.mk

LIB_CFLAGS ?= -I../../../../../../soft/common/drivers/common/utils #eigen

# LIB_CFLAGS += -I../../../../../../soft/common/drivers/linux/include #libesp, esp

# LIB_CFLAGS += -I../../../../../../soft/common/drivers/common/include #fixed_point, esp_accelerator

# LIB_CFLAGS += -I../include #gemm_stratus.h



libcrungemm.so:
	$(CXX) $(LIB_CFLAGS) c_run.cpp -shared -fPIC -o $@ $(CFLAGS) -lesp

#$(CXX) $(LIB_CFLAGS) c_run.cpp -L. $(CFLAGS) -shared -fPIC -o $@ $(LDFLAGS)

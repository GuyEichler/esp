# Copyright (c) 2011-2023 Columbia University, System Level Design Group
# SPDX-License-Identifier: Apache-2.0

EXTRA_OBJS := libcrungemm.so
CFLAGS += -L.
LIB_LDFLAGS += -L/home/geichler/Desktop/esp_accelerators/main/esp/socs/xilinx-vcu118-xcvu9p/soft-build/ariane/drivers/libesp
LDFLAGS += -lcrungemm
CXX := riscv64-unknown-linux-gnu-g++
EXTRA_HEADERS := -I/home/geichler/Desktop/esp_accelerators/main/esp/soft/common/drivers/linux/include/

APPNAME := eigen_gemm
include $(DRIVERS)/common.mk

LIB_CFLAGS ?= -I../../../../../../soft/common/drivers/common/utils #eigen

LIB_CFLAGS += -I../../../../../../soft/common/drivers/linux/include #libesp, esp

LIB_CFLAGS += -I../../../../../../soft/common/drivers/common/include #fixed_point, esp_accelerator

LIB_CFLAGS += -I../include #gemm_stratus.h



libcrungemm.so:
	$(CXX) $(LIB_CFLAGS) -c -fPIC c_run.cpp
	$(CXX) -shared -o $@ c_run.o
	$(CXX) -c -o $@ c_run.o $(CFLAGS) -lesp -lpthread -ltest -lcontig -lutils -lmonitors $(EXTRA_HEADERS)


# $(CXX) $(LIB_CFLAGS) c_run.cpp -shared -fPIC -o $@ $(CFLAGS) $(LIB_LDFLAGS) -lesp

#$(CXX) $(LIB_CFLAGS) c_run.cpp -L. $(CFLAGS) -shared -fPIC -o $@ $(LDFLAGS)

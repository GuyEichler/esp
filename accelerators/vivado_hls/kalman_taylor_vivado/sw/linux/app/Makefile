# Copyright (c) 2011-2023 Columbia University, System Level Design Group
# SPDX-License-Identifier: Apache-2.0
EXTRA_CFLAGS ?=
#APPNAME := kalman

#EXTRA_OBJS := libkalmansw.so
CFLAGS += -L.
#LDFLAGS += -lkalmansw
CXX := riscv64-unknown-linux-gnu-g++

CFLAGS += -I../../../../../common/

include $(DRIVERS)/common.mk

LIB_CFLAGS ?= -I../../../../../../soft/common/drivers/common/utils

# libkalmansw.so:
# 	$(CXX) $(LIB_CFLAGS) c_run.cpp -L. -shared -fPIC -o $@

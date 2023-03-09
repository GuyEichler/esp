// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef INC_ESPACC_CONFIG_H
#define INC_ESPACC_CONFIG_H

// User defined constants

// Data type

#define IS_TYPE_FIXED_POINT 1
#define FRAC_BITS 20
#define IS_TYPE_UINT 0
#define IS_TYPE_INT 0
#define IS_TYPE_FLOAT 0

// In/out arrays

#define SIZE_IN_CHUNK_DATA 1024

#define SIZE_OUT_CHUNK_DATA 1024

#define SIZE_OUT_BIT_DATA 1024 / DATA_BITWIDTH

#define DATA_BITWIDTH_LOG 5 //log2(32)

#define H_MAX 16
#define D_MAX 16

#endif

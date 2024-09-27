// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef INC_ESPACC_CONFIG_H
#define INC_ESPACC_CONFIG_H

// User defined constants

// Data type

#define IS_TYPE_FIXED_POINT 1
#define FRAC_BITS 32
#define WIDTH_BITS 64
#define IS_TYPE_UINT 0
#define IS_TYPE_INT 0
#define IS_TYPE_FLOAT 0

// In/out arrays

#define X_MAX 6
#define Z_MAX 164
#define CHUNK_MAX 1

#define SIZE_IN_CHUNK_DATA Z_MAX+X_MAX+X_MAX*X_MAX*3+Z_MAX*Z_MAX+Z_MAX*X_MAX

#define SIZE_OUT_CHUNK_DATA (X_MAX+X_MAX*X_MAX)*CHUNK_MAX //65792

#endif

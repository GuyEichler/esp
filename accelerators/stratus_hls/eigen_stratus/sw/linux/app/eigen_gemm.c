// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0

/* #include "libesp.h" */
#include "cfg.h"
#include "c_run.h"

int main(int argc, char **argv)
{
    unsigned in_len;
    unsigned in1_len;
    unsigned out_len;
    unsigned in_size;
    unsigned out_size;
    unsigned size;

    token_t *acc_buf;

    acc_buf           = (token_t *)esp_alloc(MAX_SIZE);
    cfg_000[0].hw_buf = acc_buf;

    unsigned *do_relu    = &(gemm_cfg_000[0].do_relu);   // filters negative features
    unsigned *transpose  = &(gemm_cfg_000[0].transpose); // If one of the matrices is transposed
    unsigned *ninputs    = &(gemm_cfg_000[0].ninputs);   // how many batches of matrix couples
    unsigned *d3         = &(gemm_cfg_000[0].d3);        // cols in rhs matrix
    unsigned *d2         = &(gemm_cfg_000[0].d2);        // cols in lhs matrix, rows in rhs matrix
    unsigned *d1         = &(gemm_cfg_000[0].d1);        // rows in lhs matrix
    unsigned *st_offset  = &(gemm_cfg_000[0].st_offset);
    unsigned *ld_offset1 = &(gemm_cfg_000[0].ld_offset1);
    unsigned *ld_offset2 = &(gemm_cfg_000[0].ld_offset2);
    unsigned *src_offset = &(gemm_cfg_000[0].src_offset);
    unsigned *dst_offset = &(gemm_cfg_000[0].dst_offset);

    int i, j, k;
    for (i = 2; i < 10; i++) {
        j = i;
        k = i;
        c_run_gemm(i, j, k, (void *)cfg_000, do_relu, transpose, ninputs, d1, d2, d3, st_offset, ld_offset1, ld_offset2,
                   src_offset, dst_offset, acc_buf);
    }

    for (i = 10; i < 100; i+=10) {
        j = i;
        k = i;
        c_run_gemm(i, j, k, (void *)cfg_000, do_relu, transpose, ninputs, d1, d2, d3, st_offset, ld_offset1, ld_offset2,
                   src_offset, dst_offset, acc_buf);
    }

    for (i = 100; i < 1000; i+=100) {
        j = i;
        k = i;
        c_run_gemm(i, j, k, (void *)cfg_000, do_relu, transpose, ninputs, d1, d2, d3, st_offset, ld_offset1, ld_offset2,
                   src_offset, dst_offset, acc_buf);
    }

    for (i = 1000; i < 10000; i+=1000) {
        j = i;
        k = i;
        c_run_gemm(i, j, k, (void *)cfg_000, do_relu, transpose, ninputs, d1, d2, d3, st_offset, ld_offset1, ld_offset2,
                   src_offset, dst_offset, acc_buf);
    }

    for (i = 10000; i < 100000; i+=10000) {
        j = i;
        k = i;
        c_run_gemm(i, j, k, (void *)cfg_000, do_relu, transpose, ninputs, d1, d2, d3, st_offset, ld_offset1, ld_offset2,
                   src_offset, dst_offset, acc_buf);
    }


    // c_run_ekf(argc, argv, (void *)cfg_000, do_relu, transpose, ninputs, d1, d2, d3, st_offset, ld_offset1,
    // ld_offset2, src_offset,
    //            dst_offset, acc_buf);

    // free
    esp_free(acc_buf);

    return 0;
}

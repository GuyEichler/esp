/* <<--functions-->> */
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* template<typename MatrixType> */
/* MatrixType operator*(const MatrixType &A, const MatrixType &B); */
void c_run_gemm(int m, int n, int p, void *x, unsigned *do_relu, unsigned *transpose, unsigned *ninputs, unsigned *d1, unsigned *d2,
                unsigned *d3, unsigned *st_offset, unsigned *ld_offset1, unsigned *ld_offset2, unsigned *src_offset,
                unsigned *dst_offset, int *acc_buf);

void c_run_ekf(int argc, char **argv, void *x, unsigned *do_relu, unsigned *transpose, unsigned *ninputs, unsigned *d1, unsigned *d2,
               unsigned *d3, unsigned *st_offset, unsigned *ld_offset1, unsigned *ld_offset2, unsigned *src_offset,
               unsigned *dst_offset, int *acc_buf);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

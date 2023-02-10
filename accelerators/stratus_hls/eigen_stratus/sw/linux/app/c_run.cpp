#include "c_run.h"
#include "eigen/Eigen/Dense"
//#include "eigen/Eigen/SVD"
#include "gemm_stratus.h"

//#include "cfg.h"
// #include <libesp.h>

#include <time.h>
#include <bits/stdc++.h>
#include <stdio.h>
#include <iostream>

/* <<--params-def-->> */
#define DO_RELU 0
#define TRANSPOSE 1
#define NINPUTS 2
#define D3 8
#define D2 8
#define D1 8
#define ST_OFFSET (NINPUTS * (D1 * D2 + D2 * D3))
#define LD_OFFSET1 0
#define LD_OFFSET2 (NINPUTS * (D1 * D2))

extern "C" void esp_dummy(void* x);


// struct gemm_stratus_access gemm_cfg_000[] = {
// 	{
// 		/* <<--descriptor-->> */
// 		.do_relu = DO_RELU,
// 		.transpose = TRANSPOSE,
// 		.ninputs = NINPUTS,
// 		.d3 = D3,
// 		.d2 = D2,
// 		.d1 = D1,
// 		.st_offset = ST_OFFSET,
// 		.ld_offset1 = LD_OFFSET1,
// 		.ld_offset2 = LD_OFFSET2,
// 		.src_offset = 0,
// 		.dst_offset = 0,
// 		.esp.coherence = ACC_COH_NONE,
// 		.esp.p2p_store = 0,
// 		.esp.p2p_nsrcs = 0,
// 		.esp.p2p_srcs = {"", "", "", ""},
// 	}
// };

// extern "C" {
// extern esp_thread_info_t cfg_000 = {
	
// 		.run = true,
// 		.devname = "gemm_stratus.0",
// 		.ioctl_req = GEMM_STRATUS_IOC_ACCESS
// 		//.esp_desc = &(gemm_cfg_000[0].esp),
	
// };
// }

using namespace Eigen;
using namespace std;

void* cfg;
unsigned* do_relu_i;
unsigned* transpose_i;
unsigned* ninputs_i;
unsigned* d3_i;
unsigned* d2_i;
unsigned* d1_i;
unsigned* st_offset_i;
unsigned* ld_offset1_i;
unsigned* ld_offset2_i;
unsigned* src_offset_i;
unsigned* dst_offset_i;
int* acc_buf_i;

#define FX_IL 16

static inline float fixed32_to_float(int value, int n_int_bits)
{
	unsigned shift_int = 0x3f800000 - 0x800000 * (32 - n_int_bits);
	float *shift = (float *) (&shift_int);

	return (*shift) * (float)value;
}

static inline int float_to_fixed32(float value, int n_int_bits)
{
	unsigned shift_int = 0x3f800000 + 0x800000 * (32 - n_int_bits);
	float *shift = (float *) &shift_int;

	return (int)(value * (*shift));
}

    template<typename MatrixType>
    MatrixType operator*(const MatrixType &A, const MatrixType &B)
    {

        // cout << "CUSTOM" << endl;

        // esp_run(cfg_000, 1);

        int rowsA = A.rows();
        int colsA = A.cols();
        int rowsB = B.rows();
        int colsB = B.cols();


        // Check if the matrices can be multiplied
        if (colsA != rowsB)
        {
            std::cerr << "Error: The number of columns in A must match the number of rows in B." << std::endl;
            return MatrixType();
        }

        *do_relu_i = 0;
        *transpose_i = 0;
        *ninputs_i = 1;

        *d1_i = rowsA;
        *d2_i = colsA;
        *d3_i = colsB;

        *ld_offset1_i = 0;
        *ld_offset2_i = 1 * rowsA * colsA;
        *st_offset_i = rowsA * colsA + rowsB * colsB;

        *src_offset_i = 0;
        *dst_offset_i = 0;

        // Allocate memory for the arrays
        int *A_data_int = &(acc_buf_i[0]);//new int[rowsA * colsA];
        int *B_data_int = &(acc_buf_i[(*ld_offset2_i)]);//new int[rowsB * colsB];
        int *C_data_int = &(acc_buf_i[(*st_offset_i)]);//new int[rowsA * colsB];

        // Copy the data from the matrices into the arrays
        for (int i = 0; i < rowsA; ++i)
            for (int j = 0; j < colsA; ++j){
                A_data_int[i * colsA + j] = float_to_fixed32(A(i, j), FX_IL);
            }
        for (int i = 0; i < rowsB; ++i)
            for (int j = 0; j < colsB; ++j){
                B_data_int[i * colsB + j] = float_to_fixed32(B(i, j), FX_IL);
            }

        // acc_buf_i = acc_buf;

        //Run gemm accelerator
        esp_dummy(cfg);

        // Copy the result from the array into a matrix
        MatrixType C_acc(rowsA, colsB);
        for (int i = 0; i < rowsA; ++i)
            for (int j = 0; j < colsB; ++j){
                C_acc(i, j) = fixed32_to_float(C_data_int[i * colsB + j], FX_IL);
            }

        return C_acc;
    }

    template<typename MatrixType>
    MatrixType CustomProductC(const MatrixType &A, const MatrixType &B)
    {

        // cout << "CUSTOM" << endl;

        // esp_run(cfg_000, 1);

        int rowsA = A.rows();
        int colsA = A.cols();
        int rowsB = B.rows();
        int colsB = B.cols();


        // Check if the matrices can be multiplied
        if (colsA != rowsB)
        {
            std::cerr << "Error: The number of columns in A must match the number of rows in B." << std::endl;
            return MatrixType();
        }


        // Allocate memory for the arrays
        double *A_data = new double[rowsA * colsA];
        double *B_data = new double[rowsB * colsB];
        double *C_data = new double[rowsA * colsB];

        // Copy the data from the matrices into the arrays
        for (int i = 0; i < rowsA; ++i)
            for (int j = 0; j < colsA; ++j){
                A_data[i * colsA + j] = A(i, j);
            }
        for (int i = 0; i < rowsB; ++i)
            for (int j = 0; j < colsB; ++j){
                B_data[i * colsB + j] = B(i, j);
            }


        // Perform the matrix multiplication using a custom algorithm
        for (int i = 0; i < rowsA; ++i)
            for (int j = 0; j < colsB; ++j)
            {
                double sum = 0.0;
                for (int k = 0; k < colsA; ++k)
                    sum += A_data[i * colsA + k] * B_data[k * colsB + j];
                C_data[i * colsB + j] = sum;
            }

        // Copy the result from the array into a matrix
        MatrixType C(rowsA, colsB);
        for (int i = 0; i < rowsA; ++i)
            for (int j = 0; j < colsB; ++j){
                C(i, j) = C_data[i * colsB + j];
            }

        // Deallocate memory for the arrays
        delete[] A_data;
        delete[] B_data;
        delete[] C_data;

        return C;
    }

// extern "C" {

    void c_run_gemm(void* x, unsigned* do_relu, unsigned* transpose,
		unsigned* ninputs, unsigned* d1, unsigned* d2, unsigned* d3,
		unsigned* st_offset, unsigned* ld_offset1, unsigned* ld_offset2,
                    unsigned* src_offset, unsigned* dst_offset, int* acc_buf)
    {

        //set global pointers to gemm accelerator configuration
        cfg = x;

        do_relu_i = do_relu;
        transpose_i = transpose;
        ninputs_i = ninputs;
        d3_i = d3;
        d2_i = d2;
        d1_i = d1;
        st_offset_i = st_offset;
        ld_offset1_i = ld_offset1;
        ld_offset2_i = ld_offset2;
        src_offset_i = src_offset;
        dst_offset_i = dst_offset;
        acc_buf_i = acc_buf;

        // MatrixXf A(2, 2);
        // A << 1, 2,
        //     3, 4;

        // MatrixXf B(2, 2);
        // B << 5, 6,
        //     7, 8;

        MatrixXf A = MatrixXf::Random(100,500);
        MatrixXf B = MatrixXf::Random(500,100);

        cout << setprecision(15) << endl;

        struct timespec startn, endn;
        struct timespec startn1, endn1;
        // time_t start, end;
        struct timespec startn2, endn2;
        // time_t start2, end2;

        cout << "Using the custom implementation and overloaded operator (gemm accelerator): " << endl;
        clock_gettime(CLOCK_MONOTONIC, &startn);
        // time(&start);

        MatrixXf C = A * B;

        clock_gettime(CLOCK_MONOTONIC, &endn);

        double time_taken = (endn.tv_sec - startn.tv_sec) * 1e9;
        time_taken = (time_taken + (endn.tv_nsec - startn.tv_nsec)) * 1e-9;
        // float time_taken = double(end - start);

        cout << fixed << "Time taken: " << time_taken << " sec\n" << endl;

        cout << "Using the custom implementation in regular C: " << endl;
        clock_gettime(CLOCK_MONOTONIC, &startn1);
        // time(&start);

        MatrixXf C1 = CustomProductC(A, B);

        clock_gettime(CLOCK_MONOTONIC, &endn1);
        // time(&end);

        // cout << "The product of A and B is:\n" << C << endl;

        double time_taken1 = (endn1.tv_sec - startn1.tv_sec) * 1e9;
        time_taken1 = (time_taken1 + (endn1.tv_nsec - startn1.tv_nsec)) * 1e-9;
        // float time_taken = double(end - start);

        cout << fixed << "Time taken: " << time_taken1 << " sec\n" << endl;
        // cout << endn.tv_sec << endl;
        // cout << startn.tv_sec << endl;

        cout << "Using the Eigen implementation and original operator: " << endl;
        clock_gettime(CLOCK_MONOTONIC, &startn2);
        // time(&start2);

        MatrixXf C2 = A.operator*(B);

        clock_gettime(CLOCK_MONOTONIC, &endn2);
        // time(&end2);


        double time_taken2 = (endn2.tv_sec - startn2.tv_sec) * 1e9;
        time_taken2 = (time_taken2 + (endn2.tv_nsec - startn2.tv_nsec)) * 1e-9;
        // float time_taken2 = double(end2 - start2);

        // cout << "The product of A and B is:\n" << C2 << endl;

        cout << fixed << "Time taken: " << time_taken2 << " sec" << endl;
        // cout << endn2.tv_sec << endl;
        // cout << startn2.tv_sec << endl;

        // Eigen::MatrixXf C2 = Eigen::MatrixBase<Eigen::MatrixXf>::operator*(A,B);

        //return 0;

        MatrixXf C01 = C - C1;

        cout << endl;
        cout << "Maximum difference between accelerator and C: " << C01.maxCoeff() << endl;
        cout << "Minimum difference between accelerator and C: " << C01.minCoeff() << endl;

    }

    // void printR()
    // {
    //     cout <<"\n+ Software matrix R is:\n" << R << endl;
    // }

// } /* extern "C" */

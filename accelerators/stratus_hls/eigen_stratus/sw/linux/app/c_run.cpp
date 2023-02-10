#include "c_run.h"
#include "eigen/Eigen/Dense"
//#include "eigen/Eigen/SVD"
#include "gemm_stratus.h"

//#include "cfg.h"
//#include <libesp.h>

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

// extern "C" void esp_dummy();

using namespace Eigen;
using namespace std;

    template<typename MatrixType>
    MatrixType operator*(const MatrixType &A, const MatrixType &B)
    {

        // cout << "CUSTOM" << endl;

        // esp_dummy();

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
            for (int j = 0; j < colsA; ++j)
                A_data[i * colsA + j] = A(i, j);
        for (int i = 0; i < rowsB; ++i)
            for (int j = 0; j < colsB; ++j)
                B_data[i * colsB + j] = B(i, j);

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
            for (int j = 0; j < colsB; ++j)
                C(i, j) = C_data[i * colsB + j];

        // Deallocate memory for the arrays
        delete[] A_data;
        delete[] B_data;
        delete[] C_data;

        return C;
    }

// extern "C" {

    void c_run_gemm()
    {
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
        // time_t start, end;
        struct timespec startn2, endn2;
        // time_t start2, end2;

        cout << "Using the custom implementation and overloaded operator: " << endl;
        clock_gettime(CLOCK_MONOTONIC, &startn);
        // time(&start);

        MatrixXf C = A * B;

        clock_gettime(CLOCK_MONOTONIC, &endn);
        // time(&end);

        // cout << "The product of A and B is:\n" << C << endl;

        double time_taken = (endn.tv_sec - startn.tv_sec) * 1e9;
        time_taken = (time_taken + (endn.tv_nsec - startn.tv_nsec)) * 1e-9;
        // float time_taken = double(end - start);

        cout << fixed << "Time taken: " << time_taken  << " sec\n" << endl;
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

    }

    // void printR()
    // {
    //     cout <<"\n+ Software matrix R is:\n" << R << endl;
    // }

// } /* extern "C" */

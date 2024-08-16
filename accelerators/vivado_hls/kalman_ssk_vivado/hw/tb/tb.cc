// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "../inc/espacc_config.h"
#include "../inc/espacc.h"
#include "../inc/espacc_tb_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <string>

#include "A_array.h"
#include "H_array.h"
#include "W_array.h"
#include "Q_array.h"
#include "initial_state_array.h"
#include "measurements_array.h"
#include "prediction_array.h"
#include "real_array.h"
#include "P_array.h"

// #include "A_array_soma.h"
// #include "H_array_soma.h"
// #include "W_array_soma.h"
// #include "Q_array_soma.h"
// #include "initial_state_array_soma.h"
// #include "measurements_array_soma.h"
// #include "prediction_array_soma.h"
// #include "real_array_soma.h"
// #include "P_array_soma.h"

// #include "A_array_hc.h"
// #include "H_array_hc.h"
// #include "W_array_hc.h"
// #include "Q_array_hc.h"
// #include "initial_state_array_hc.h"
// #include "measurements_array_hc.h"
// #include "prediction_array_hc.h"
// #include "real_array_hc.h"
// #include "P_array_hc.h"

#include <eigen3/Eigen/Dense>
using namespace Eigen;

#define STATES 6
#define NEURONS 164
#define TIME_STAMPS 100
#define CHUNKS 1
#define BATCHES TIME_STAMPS / CHUNKS

// #define STATES 6
// #define NEURONS 52
// #define TIME_STAMPS 500
// #define CHUNKS 100
// #define BATCHES TIME_STAMPS / CHUNKS

// #define STATES 6
// #define NEURONS 46
// #define TIME_STAMPS 100
// #define CHUNKS 100
// #define BATCHES TIME_STAMPS / CHUNKS

#define DOWNSMPL 1
#define NOISE 0
#define SCALE 1

#define SSK_TIME_STAMPS 400

int main(int argc, char **argv) {

    printf("****start*****\n");

    /* <<--params-->> */
    const unsigned inv_reset = 1;
    const unsigned inv_num = 2;
    const unsigned chunks = CHUNKS;
    const unsigned iter = BATCHES;
    const unsigned x_dim = STATES;
    const unsigned z_dim = NEURONS;

    uint32_t in_words_adj;
    uint32_t in_words_adj_z;
    uint32_t out_words_adj;
    uint32_t in_size;
    uint32_t out_size;
    uint32_t dma_in_size;
    uint32_t dma_out_size;
    uint32_t dma_size;

    //time_t t;
    //srand((unsigned) time(&t));
    srand(123);

    in_words_adj = round_up(x_dim + x_dim * x_dim * 3 + 2 * z_dim * z_dim + z_dim * x_dim, VALUES_PER_WORD);
    in_words_adj_z = round_up(z_dim * chunks, VALUES_PER_WORD);
    out_words_adj = round_up((x_dim + x_dim * x_dim) * chunks, VALUES_PER_WORD);
    in_size = in_words_adj + in_words_adj_z * iter;
    out_size = out_words_adj * (iter);

    dma_in_size = in_size / VALUES_PER_WORD;
    dma_out_size = out_size / VALUES_PER_WORD;
    dma_size = dma_in_size + dma_out_size;

    dma_word_t *mem=(dma_word_t*) malloc(dma_size * sizeof(dma_word_t));
    word_t *inbuff=(word_t*) malloc(in_size * sizeof(word_t));
    word_t *outbuff=(word_t*) malloc(out_size * sizeof(word_t));
    word_t *outbuff_gold= (word_t*) malloc(out_size * sizeof(word_t));
    dma_info_t load;
    dma_info_t store;

    // Prepare input data

    //Compute golden P with Eigen
    Matrix<float, NEURONS, STATES> Mat_H = Map<Matrix<float, NEURONS, STATES, RowMajor> >(H);
    Matrix<float, STATES, STATES> Mat_F = Map<Matrix<float, STATES, STATES, RowMajor> >(A);
    Matrix<float, STATES, STATES> Mat_P = Matrix<float, STATES, STATES>::Zero();
    Matrix<float, STATES, STATES> Mat_Q = Map<Matrix<float, STATES, STATES, RowMajor> >(W);
    Matrix<float, NEURONS, NEURONS> Mat_R = Map<Matrix<float, NEURONS, NEURONS, RowMajor> >(Q);
    Matrix<float, STATES, STATES> Mat_I = Matrix<float, STATES, STATES>::Identity();
    Matrix<float, NEURONS, NEURONS> prev_S_inv;
    Matrix<float, NEURONS, NEURONS> prev_S_inv_iter;
    Matrix<float, STATES, STATES> Mat_P_iter = Matrix<float, STATES, STATES>::Zero();
    Matrix<float, STATES, STATES> Mat_P_gauss = Matrix<float, STATES, STATES>::Zero();

    //prepare SSK_inv
    Matrix<float, NEURONS, NEURONS> Mat_SSK_inv;
    float SSK_inv[Z_MAX * Z_MAX];

    for(int i = 0; i < SSK_TIME_STAMPS; i++)
    {
        Matrix<float, NEURONS, NEURONS> Mat_S = Mat_H * (Mat_F * Mat_P * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() + Mat_R;
        Matrix<float, STATES, NEURONS> Mat_K = (Mat_F * Mat_P * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() * Mat_S.inverse();
        Mat_P = (Mat_I -  Mat_K * Mat_H) * (Mat_F * Mat_P * Mat_F.transpose() + Mat_Q);
        // if(i == 0)
        //     std::cout << "Matrix P = " << Mat_P << std::endl;

        if(i < TIME_STAMPS)
        {
            float P_flat[STATES*STATES];
            Map<Matrix<float, STATES, STATES, RowMajor> >(&P_flat[0], Mat_P.rows(), Mat_P.cols()) = Mat_P;

            for(unsigned j = 0; j < x_dim + x_dim * x_dim; j++)
            {
                if(j < x_dim)
                {
                    outbuff_gold[i/DOWNSMPL * out_words_adj/chunks + j] = (word_t) prediction[STATES * (i+1) + j];
                }
                else
                {
                    outbuff_gold[i/DOWNSMPL * out_words_adj/chunks + j] = (word_t) P_flat[j - x_dim];
                }
            }
        }

        //Last iteration save SSK
        if(i == SSK_TIME_STAMPS - 1)
            Mat_SSK_inv = Mat_S.inverse();

    }

    for(int i = 0; i < NEURONS; i++)
        for(int j = 0; j < NEURONS; j++)
            SSK_inv[i * NEURONS + j] = Mat_SSK_inv(i, j);

    // for(unsigned i = 0; i < iter; i++)
    //     for(unsigned j = 0; j < z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim; j++)
    //         inbuff[i * in_words_adj + j] = (word_t) j;

    for(unsigned i = 0; i < iter; i++)
    {//z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim
        unsigned j = 0;

        // //Z
        // for(; j < z_dim; j++)
        // {
        //     if(i == 0)
        //         inbuff[i * in_words_adj + j] = (word_t) measurements[NEURONS * (i+1) + j];
        //     else
        //         inbuff[in_words_adj + (i-1) * in_words_adj_z + j] = (word_t) measurements[NEURONS * (i+1) + j];
        //     // if(i == 1)
        //     //      printf("Value of Z = %f index %d\n", measurements[NEURONS * (i+1) + j], in_words_adj + (i-1) * in_words_adj_z + j);
        // }

        if(i == 0) //only for first iteration
        {
            //X
            for(; j < x_dim; j++)
            {
                inbuff[i * in_words_adj + j] = (word_t) initial[j];
            }

            //P
            for(; j < x_dim + x_dim * x_dim; j++)
            {
                inbuff[i * in_words_adj + j] = (word_t) 0.0;
            }

            //F
            for(; j < x_dim + x_dim * x_dim * 2; j++)
            {

                float rand_val = (float)rand() / RAND_MAX;
                rand_val = 2.0 * NOISE * rand_val - NOISE;

                inbuff[i * in_words_adj + j] = (word_t) A[j - (x_dim + x_dim * x_dim)];// + rand_val;
            }

            //Q
            for(; j < x_dim + x_dim * x_dim * 3; j++)
            {

                float rand_val = (float)rand() / RAND_MAX;
                rand_val = 2.0 * NOISE * rand_val - NOISE;

                inbuff[i * in_words_adj + j] = (word_t) W[j - (x_dim + x_dim * x_dim * 2)];// + rand_val;
            }

            //R
            for(; j < x_dim + x_dim * x_dim * 3 + z_dim * z_dim; j++)
            {

                float rand_val = (float)rand() / RAND_MAX;
                rand_val = 2.0 * NOISE * rand_val - NOISE;

                inbuff[i * in_words_adj + j] = (word_t) Q[j - (x_dim + x_dim * x_dim * 3)];// + rand_val;
            }

            //H
            for(; j < x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim; j++)
            {

                float rand_val = (float)rand() / RAND_MAX;
                rand_val = 2.0 * NOISE * rand_val - NOISE;

                inbuff[i * in_words_adj + j] = (word_t) H[j - (x_dim + x_dim * x_dim * 3 + z_dim * z_dim)];// + rand_val;
                //printf("Value of H = %f \n", measurements[NEURONS * (i+1) + j]);
            }

            //SSK
            for(; j < x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim + z_dim * z_dim; j++)
            {
                inbuff[i * in_words_adj + j] = (word_t) SSK_inv[j - (x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim)];// + rand_val;
                //printf("Value of H = %f \n", measurements[NEURONS * (i+1) + j]);
            }
        }

        unsigned base_index = (x_dim + x_dim * x_dim * 3 + 2 * z_dim * z_dim + z_dim * x_dim);
        //Z

        unsigned j2 = j;

        if (i == 0)
            for(; j < base_index + z_dim * chunks; j++)
            {

                float rand_val = (float)rand() / RAND_MAX;
                rand_val = 2.0 * NOISE * rand_val - NOISE;

                inbuff[j] = (word_t) measurements[NEURONS * (i+1) + j2 - base_index] * SCALE + rand_val;

                //printf("Value of j = %d, j2 = %d , base = %d, mod = %d \n", j, j2, base_index, (j-base_index)%z_dim);

                if((j+1 - base_index) != 0 && ((j+1 - base_index) % z_dim == 0) && DOWNSMPL > 1)
                {
                    j2 += z_dim * (DOWNSMPL - 1) + 1;
                }
                else
                    j2++;

                // printf("Value of Z = %f index %d inbuff %d\n", measurements[NEURONS * (i+1) + j - base_index], NEURONS * (i+1) + j - base_index, j);
            }
        else
            for(; j < z_dim * chunks; j++)
            {

                float rand_val = (float)rand() / RAND_MAX;
                rand_val = 2.0 * NOISE * rand_val - NOISE;

                inbuff[in_words_adj + i * in_words_adj_z + j] = (word_t) measurements[NEURONS * (i*chunks+1) + j2] * SCALE + rand_val;

                if(j+1 != 0 && (j+1  %  z_dim == 0) && DOWNSMPL > 1)
                {
                    j2 += z_dim * (DOWNSMPL - 1) + 1;
                }
                else
                    j2++;

                // printf("Value of Z = %f index %d inbuff %d i = %d\n", measurements[NEURONS * (i*chunks+1) + j], NEURONS * (i*chunks+1) + j, in_words_adj + i * in_words_adj_z + j ,i);
            }
    }

    for(unsigned i = 0; i < dma_in_size; i++)
	for(unsigned k = 0; k < VALUES_PER_WORD; k++)
	    mem[i].word[k] = inbuff[i * VALUES_PER_WORD + k];

    // Set golden output
    // for(unsigned i = 0; i < iter; i++)
    //     for(unsigned j = 0; j < x_dim + x_dim * x_dim; j++)
    //         outbuff_gold[i * out_words_adj + j] = (word_t) j;


    //Files for printing values
    std::ofstream outfile_eigen("norm_eigen.h");
    std::ofstream outfile_gauss("norm_gauss.h");
    std::string norm_newton = "norm_newton" + std::to_string(inv_num) + std::to_string(inv_reset) + ".h";
    std::ofstream outfile_newton(norm_newton);
    std::string norm_newton_final = "norm_newton_final" + std::to_string(inv_num) + std::to_string(inv_reset) + ".h";
    std::ofstream outfile_newton_final(norm_newton_final);

    std::ofstream outfile_originalX("original_positionX.h");
    std::ofstream outfile_originalY("original_positionY.h");
    std::ofstream outfile_sotaX("sota_positionX.h");
    std::ofstream outfile_sotaY("sota_positionY.h");
    std::string predictionX = "predicted_positionX" + std::to_string(inv_num) + std::to_string(inv_reset) + ".h";
    std::string predictionY = "predicted_positionY" + std::to_string(inv_num) + std::to_string(inv_reset) + ".h";
    std::ofstream outfile_predictionX(predictionX);
    std::ofstream outfile_predictionY(predictionY);

    int counter = 1;

    // if(outfile_eigen.is_open()){
    //     outfile_eigen << "norm_eigen = [";
    // }

    // if(outfile_gauss.is_open()){
    //     outfile_gauss << "norm_gauss = [";
    // }

    // if(outfile_newton.is_open()){
    //     outfile_newton << "norm_newton = [";
    // }

    // for(int i = 0; i < TIME_STAMPS * DOWNSMPL; i+=DOWNSMPL)
    // {
    //     Matrix<float, NEURONS, NEURONS> Mat_S = Mat_H * (Mat_F * Mat_P * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() + Mat_R;
    //     Matrix<float, STATES, NEURONS> Mat_K = (Mat_F * Mat_P * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() * Mat_S.inverse();
    //     Mat_P = (Mat_I -  Mat_K * Mat_H) * (Mat_F * Mat_P * Mat_F.transpose() + Mat_Q);
    //     // if(i == 0)
    //     //     std::cout << "Matrix P = " << Mat_P << std::endl;

    //     float P_flat[STATES*STATES];
    //     Map<Matrix<float, STATES, STATES, RowMajor> >(&P_flat[0], Mat_P.rows(), Mat_P.cols()) = Mat_P;

    //     for(unsigned j = 0; j < x_dim + x_dim * x_dim; j++)
    //     {
    //         if(j < x_dim)
    //         {
    //             outbuff_gold[i/DOWNSMPL * out_words_adj/chunks + j] = (word_t) prediction[STATES * (i+1) + j];
    //             // if(j == 2){
    //             //     if(outfile_originalX.is_open()){
    //             //         float x = real_out[STATES * (i+1) + j];
    //             //         outfile_originalX << x << std::endl;
    //             //         float x_pred = prediction[STATES * (i+1) + j];
    //             //         outfile_sotaX << x_pred << std::endl;
    //             //     }
    //             // }
    //             // if(j == 3){
    //             //     if(outfile_originalY.is_open()){
    //             //         float y = real_out[STATES * (i+1) + j];
    //             //         outfile_originalY << y << std::endl;
    //             //         float y_pred = prediction[STATES * (i+1) + j];
    //             //         outfile_sotaY << y_pred << std::endl;
    //             //     }
    //             // }
    //         }
    //         else
    //         {
    //             outbuff_gold[i/DOWNSMPL * out_words_adj/chunks + j] = (word_t) P_flat[j - x_dim];
    //         }
    //     }

    //     Matrix<float, NEURONS, NEURONS> Mat_S_iter;

        // //Norm check

        // //I - A*X_0
        // Matrix<float, NEURONS, NEURONS> res = Mat_S * Mat_S.inverse();
        // Matrix<float, NEURONS, NEURONS> I_res = Matrix<float, NEURONS, NEURONS>::Identity();
        // Matrix<float, NEURONS, NEURONS> final_res = I_res - res;
        // float norm = final_res.norm();
        // std::cout << "Iteration: " << i << " Norm for calculation = " << norm << std::endl;

        // if(outfile_eigen.is_open()){
        //     outfile_eigen << norm << std::endl;
        //     // if(i != TIME_STAMPS - 1){
        //     //     outfile_eigen << ", ";
        //     // }
        //     // else{
        //     //     outfile_eigen << "]";
        //     //     outfile_eigen.close();
        //     // }
        // }

        // //Manual gauss inverse
        // Matrix<float, NEURONS, NEURONS> Mat_S_gauss = Mat_H * (Mat_F * Mat_P_gauss * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() + Mat_R;

        // float S_array_gauss[NEURONS][NEURONS] = {{0}};
        // float S_inv_array_gauss[NEURONS][NEURONS] = {{0}};

        // for(int k = 0; k < NEURONS; k++)
        //     for(int d = 0; d < NEURONS; d++)
        //     {
        //         //S_inv_array_gauss[k][d] = S_inv(k,d);
        //         S_array_gauss[k][d] = Mat_S_gauss(k,d);
        //         if(k == d)
        //             S_inv_array_gauss[k][d] = 1.0;
        //     }

        // int inv_ok = inverse_tb<NEURONS, float>(S_array_gauss, S_inv_array_gauss, z_dim);

        // Map<Matrix<float, NEURONS, NEURONS, RowMajor> > Mat_S_inv_gauss(&S_inv_array_gauss[0][0], NEURONS, NEURONS);

        // //std::cout << "Eigen: " << Mat_S.inverse() << std::endl;

        // // for(int k = 0; k < NEURONS; k++)
        // // {
        // //     std::cout << std::endl;
        // //     for(int d = 0; d < NEURONS; d++)
        // //     {
        // //         std::cout << S_inv_array_gauss[k][d] <<" ";
        // //     }
        // // }

        // // std::cout << std::endl;

        // //std::cout << "Gauss: " << Mat_S_inv_gauss << std::endl;

        // Matrix<float, STATES, NEURONS> Mat_K_gauss = (Mat_F * Mat_P_gauss * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() * Mat_S_inv_gauss;
        // Mat_P_gauss = (Mat_I -  Mat_K_gauss * Mat_H) * (Mat_F * Mat_P_gauss * Mat_F.transpose() + Mat_Q);


        // //Check gauss norm
        // Matrix<float, NEURONS, NEURONS> res_gauss = Mat_S_gauss * Mat_S_inv_gauss;
        // Matrix<float, NEURONS, NEURONS> I_res_gauss = Matrix<float, NEURONS, NEURONS>::Identity();
        // Matrix<float, NEURONS, NEURONS> final_res_gauss = I_res_gauss - res_gauss;
        // float norm_gauss = final_res_gauss.norm();
        // Matrix<float, NEURONS, NEURONS> diff_mat_S_g = Mat_S_gauss - Mat_S;
        // float diff_g = std::max(diff_mat_S_g.maxCoeff(), (-1)*diff_mat_S_g.minCoeff());
        // //std::cout << "Iteration: " << i << " Norm for gauss = " << norm_gauss << std::endl;
        // //std::cout << "Iteration: " << i << " max diff from gauss = " << diff_g << std::endl;
        // // JacobiSVD<MatrixXf> svd(Mat_S_gauss);
        // // VectorXf singularValues = svd.singularValues();
        // // double cond_spectral = singularValues.maxCoeff() / singularValues.minCoeff();
        // // std::cout << "Condition for gauss: " << cond_spectral << std::endl;
        // // double det_gauss = Mat_S_gauss.determinant();;
        // // std::cout << "Determinant gauss: " << det_gauss << std::endl;

        // if(outfile_gauss.is_open()){
        //     outfile_gauss << norm_gauss << std::endl;
        //     // if(i != TIME_STAMPS - 1){
        //     //     outfile_gauss << ", ";
        //     // }
        //     // else{
        //     //     outfile_gauss << "]";
        //     //     outfile_gauss.close();
        //     // }
        // }

        // //If we are in the first time stamp the norm of newton is the same as gauss
        // if(i == 0){
        //     if(outfile_newton.is_open()){
        //         outfile_newton << norm_gauss << std::endl;
        //         outfile_newton_final << norm_gauss << std::endl;
        //         // if(i != TIME_STAMPS - 1){
        //         //     outfile_newton << ", ";
        //         // }
        //         // else{
        //         //     outfile_newton << "]";
        //         //     outfile_newton.close();
        //         // }
        //     }
        // }


        // //Iterative inverse
        // if(i > 0)
        // {

        //     Mat_S_iter = Mat_H * (Mat_F * Mat_P_iter * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() + Mat_R;

        //     //Check iterative norm
        //     //Calculate with the iterative norm
        //     Matrix<float, NEURONS, NEURONS> res_iter = Mat_S_iter * prev_S_inv_iter;
        //     Matrix<float, NEURONS, NEURONS> I_res_iter = Matrix<float, NEURONS, NEURONS>::Identity();
        //     Matrix<float, NEURONS, NEURONS> final_res_iter = I_res_iter - res_iter;
        //     float norm_iter = final_res_iter.norm();
        //     // Matrix<float, NEURONS, NEURONS> diff_mat_S = Mat_S_iter - Mat_S;
        //     // float diff = std::max(diff_mat_S.maxCoeff(), (-1)*diff_mat_S.minCoeff());
        //     //std::cout << "Iteration: " << i << " Norm for approximation = " << norm_iter << std::endl;
        //     //std::cout << "Iteration: " << i << " max diff = " << diff << std::endl;

        //     if(counter < inv_reset || inv_reset == 0){
        //         if(outfile_newton.is_open()){
        //             outfile_newton << norm_iter << std::endl;
        //         }
        //     }

        //     //Initial condition
        //     // //Set the initial condition to the previous gauss
        //     // Matrix<float, NEURONS, NEURONS> S_inv =  prev_S_inv;
        //     //Set the initial condition to the previous newton
        //     Matrix<float, NEURONS, NEURONS> S_inv =  prev_S_inv_iter;

        //     //std::cout << "prev S inv: " << S_inv << std::endl;

        //     float S_inv_array[NEURONS][NEURONS] = {{0}};
        //     float Mat_S_array[NEURONS][NEURONS] = {{0}};
        //     //Map<Matrix<float, NEURONS, NEURONS, RowMajor> > iterative_out(&next_out[0][0], N, N);
        //     Matrix<float, NEURONS, NEURONS> S_inv_final;
        //     float S_inv_final_array[NEURONS][NEURONS] = {{0}};
        //     Matrix<float, NEURONS, NEURONS> S_inv_final2;
        //     float S_inv_final2_array[NEURONS][NEURONS] = {{0}};
        //     for(int k = 0; k < NEURONS; k++)
        //         for(int d = 0; d < NEURONS; d++)
        //         {
        //             S_inv_array[k][d] = S_inv(k,d);
        //             Mat_S_array[k][d] = Mat_S_iter(k,d);
        //         }

        //     //Map<Matrix<float, NEURONS, NEURONS, RowMajor> > iterative_out(&next_out[0][0], N, N);

        //     for(unsigned l = 0; l < inv_num; l++){
        //         if(l == 0){
        //             iterative_inverse_tb<NEURONS>(Mat_S_array, S_inv_array, S_inv_final_array, z_dim);
        //             //std::cout << "Iterative 0" << std::endl;
        //         }
        //         else if(l % 2 == 1){
        //             iterative_inverse_tb<NEURONS>(Mat_S_array, S_inv_final_array, S_inv_final2_array, z_dim);
        //             //std::cout << "Iterative 1" << std::endl;
        //         }
        //         else{// if(i % 2 == 0){
        //             iterative_inverse_tb<NEURONS>(Mat_S_array, S_inv_final2_array, S_inv_final_array, z_dim);
        //             //std::cout << "Iterative 2" << std::endl;
        //         }
        //     }

        //     if(counter < inv_reset || inv_reset == 0)
        //     {
        //         if(inv_num % 2 == 0)
        //         {//final2
        //             Map<Matrix<float, NEURONS, NEURONS, RowMajor> > iterative_out(&S_inv_final2_array[0][0], NEURONS, NEURONS);
        //             prev_S_inv_iter = iterative_out;
        //             //std::cout << "Set final 2" << std::endl;
        //         }
        //         else
        //         {//final
        //             Map<Matrix<float, NEURONS, NEURONS, RowMajor> > iterative_out(&S_inv_final_array[0][0], NEURONS, NEURONS);
        //             prev_S_inv_iter = iterative_out;
        //             //std::cout << "Set final 1" << std::endl;
        //         }

        //         Matrix<float, STATES, NEURONS> Mat_K_iter = (Mat_F * Mat_P_iter * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() * prev_S_inv_iter;//Mat_S.inverse();
        //         Mat_P_iter = (Mat_I -  Mat_K_iter * Mat_H) * (Mat_F * Mat_P_iter * Mat_F.transpose() + Mat_Q);

        //         //Calculate with the final iterative norm
        //         Matrix<float, NEURONS, NEURONS> res_iter_final = Mat_S_iter * prev_S_inv_iter;
        //         Matrix<float, NEURONS, NEURONS> I_res_iter_final = Matrix<float, NEURONS, NEURONS>::Identity();
        //         Matrix<float, NEURONS, NEURONS> final_res_iter_final = I_res_iter_final - res_iter_final;
        //         float norm_iter_final = final_res_iter_final.norm();

        //         // JacobiSVD<MatrixXf> svd_iter(Mat_S_iter);
        //         // VectorXf singularValues_iter = svd_iter.singularValues();
        //         // double cond_spectral_iter = singularValues_iter.maxCoeff() / singularValues_iter.minCoeff();
        //         // std::cout << "Condition for newton: " << cond_spectral_iter
        //         // << std::endl;
        //         // double det_iter = Mat_S_iter.determinant();;
        //         // std::cout << "Determinant newton: " << det_iter << std::endl;

        //         if(outfile_newton.is_open()){
        //             outfile_newton_final << norm_iter_final << std::endl;
        //         }

        //         counter++;
        //     }
        //     else
        //     {
        //         //std::cout << "HERE" << std::endl;
        //         Mat_S_iter = Mat_H * (Mat_F * Mat_P_iter * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() + Mat_R;

        //         //Manual gauss inverse
        //         float S_array_gauss_iter[NEURONS][NEURONS] = {{0}};
        //         float S_inv_array_gauss_iter[NEURONS][NEURONS] = {{0}};

        //         for(int k = 0; k < NEURONS; k++)
        //             for(int d = 0; d < NEURONS; d++)
        //             {
        //                 //S_inv_array_gauss[k][d] = S_inv(k,d);
        //                 S_array_gauss_iter[k][d] = Mat_S_iter(k,d);
        //                 if(k == d)
        //                     S_inv_array_gauss_iter[k][d] = 1.0;
        //             }

        //         int inv_ok = inverse_tb<NEURONS, float>(S_array_gauss_iter, S_inv_array_gauss_iter, z_dim);

        //         Map<Matrix<float, NEURONS, NEURONS, RowMajor> > Mat_S_inv_gauss_iter(&S_inv_array_gauss[0][0], NEURONS, NEURONS);
        //         prev_S_inv_iter = Mat_S_inv_gauss_iter;

        //         Matrix<float, STATES, NEURONS> Mat_K_iter = (Mat_F * Mat_P_iter * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() * prev_S_inv_iter;//Mat_S.inverse();
        //         Mat_P_iter = (Mat_I -  Mat_K_iter * Mat_H) * (Mat_F * Mat_P_iter * Mat_F.transpose() + Mat_Q);

        //         Matrix<float, NEURONS, NEURONS> res_iter_gauss = Mat_S_iter * prev_S_inv_iter;
        //         Matrix<float, NEURONS, NEURONS> I_res_iter_gauss = Matrix<float, NEURONS, NEURONS>::Identity();
        //         Matrix<float, NEURONS, NEURONS> final_res_iter_gauss = I_res_iter_gauss - res_iter_gauss;
        //         float norm_iter_gauss = final_res_iter_gauss.norm();

        //         // JacobiSVD<MatrixXf> svd_iter_gauss(Mat_S_iter);
        //         // VectorXf singularValues_iter_gauss = svd_iter_gauss.singularValues();
        //         // double cond_spectral_iter_gauss = singularValues_iter_gauss.maxCoeff() / singularValues_iter_gauss.minCoeff();
        //         // std::cout << "Condition for newton (gauss): " << cond_spectral_iter_gauss << std::endl;
        //         // double det_iter_gauss = Mat_S_iter.determinant();
        //         // std::cout << "Determinant newton: " << det_iter_gauss << std::endl;

        //         if(outfile_newton.is_open()){
        //             outfile_newton << norm_iter_gauss << std::endl;
        //             outfile_newton_final << norm_iter_gauss << std::endl;
        //         }

        //         counter = 1;
        //     }
        // }
        // else
        // {
        //     //prev_S_inv_iter = Mat_S.inverse();
        //     prev_S_inv_iter = Mat_S_inv_gauss;
        //     //Mat_P_iter = Mat_P;
        //     Mat_P_iter = Mat_P_gauss;
        // }


        // //Update prev_S_inv
        // prev_S_inv = Mat_S.inverse();

        // // //Update prev_S_inv
        // // prev_S_inv_gauss = Mat_S.inverse();

    // }

    // for(int i = 0; i < iter; i++)
    // for(unsigned j = 0; j < x_dim + x_dim * x_dim; j++)
    // {
    //     if(j < x_dim)
    //         outbuff_gold[i * out_words_adj + j] = (word_t) prediction[STATES * (i+1) + j];
    //     else
    //     {
    //         outbuff_gold[i * out_words_adj + j] = (word_t) P_flat[i * x_dim * x_dim + j - x_dim];
    //     }
    // }

    // for(int i = 0; i < STATES*STATES; i++)
    //     printf("P_flat[%d] = %.10f\n", i, P_flat[i]);


    // for(unsigned i = 0; i < iter; i++)
    //     for(unsigned j = 0; j < x_dim + x_dim * x_dim; j++)
    //     {
    //         if(j < x_dim)
    //             outbuff_gold[i * out_words_adj + j] = (word_t) prediction[STATES * (i+1) + j];
    //         else
    //         {
    //             outbuff_gold[i * out_words_adj + j] = (word_t) P_flat[j - x_dim];
    //         }
    //     }


    // Call the TOP function
    top(mem, mem,
        /* <<--args-->> */
                 inv_reset,
                 inv_num,
         	 chunks,
         	 iter,
         	 x_dim,
         	 z_dim,
        load, store);

    // Validate
    uint32_t out_offset = dma_in_size;
    for(unsigned i = 0; i < dma_out_size; i++)
	for(unsigned k = 0; k < VALUES_PER_WORD; k++)
	    outbuff[i * VALUES_PER_WORD + k] = mem[out_offset + i].word[k];

    int errors = 0;
    float MSE = 0;
    float MAE = 0;
    float max_diff = 0;
    float thresh = 0.3;

    for(unsigned i = 0; i < iter; i++)
        for(unsigned j = 0; j < (x_dim + x_dim * x_dim)*chunks; j++)
        {
            word_t gold_val = outbuff_gold[i * out_words_adj + j];
            word_t acc_val = outbuff[i * out_words_adj + j];
            // printf("Accelerator value: %f Golden value: %f index: %d\n", acc_val, gold_val, i * out_words_adj + j);

            word_t diff = std::abs(gold_val - acc_val);
            MAE += diff;
            MSE += std::pow(diff, 2.0);
            word_t norm_diff = diff/acc_val > diff/gold_val ? diff/acc_val : diff/gold_val;
            if(norm_diff > max_diff)
                max_diff = norm_diff;

            if(j%(x_dim + x_dim*x_dim) < x_dim){
                //std::cout << "OKAY: X Accelerator value: " << acc_val << " Golden value: " << gold_val << " index: " << i * out_words_adj + j << " Iter " << i << " diff: " << diff << std::endl;

                // if(j%(x_dim + x_dim*x_dim) == 2){
                //     if(outfile_predictionX.is_open()){
                //         word_t x = acc_val;
                //         outfile_predictionX << x << std::endl;
                //     }
                // }
                // else if(j%(x_dim + x_dim*x_dim) == 3){
                //     if(outfile_predictionY.is_open()){
                //         word_t y = acc_val;
                //         outfile_predictionY << y << std::endl;
                //     }
                // }
            }
            // else
            //     std::cout << "P Accelerator value: " << acc_val << " Golden value: " << gold_val << " index: " << i * out_words_adj + j << " Iter " << i << " diff: " << diff << std::endl;

	    if (outbuff[i * out_words_adj + j] != outbuff_gold[i * out_words_adj + j]){
                if(diff/gold_val > thresh || diff/acc_val > thresh || diff/gold_val < (-1)*thresh || diff/acc_val < (-1) * thresh){
                    //printf("Accelerator value: %f Golden value: %f index: %d\n", acc_val, gold_val, i * out_words_adj + j);
                    if(j < x_dim)
                        std::cout << "BAD: X Accelerator value: " << acc_val << " Golden value: " << gold_val << " index: " << i * out_words_adj + j << " Iter " << i << " diff: " << diff << std::endl;
                    // else
                    //     std::cout << "P Accelerator value: " << acc_val << " Golden value: " << gold_val << " index: " << i * out_words_adj + j << " Iter " << i << " diff: " << diff << std::endl;

                    errors++;
                }
            }
        }

    MSE /= ((x_dim + x_dim * x_dim) * iter*chunks);
    MAE /= ((x_dim + x_dim * x_dim) * iter*chunks);
    std::cout << "Output MSE: " << MSE << std::endl;
    std::cout << "Output MAE: " << MAE << std::endl;
    std::cout << "Output max_diff: " << max_diff << std::endl;

    if (errors)
	std::cout << "Test FAILED with " << errors << " errors." << std::endl;
    else
	std::cout << "Test PASSED." << std::endl;

    // Free memory

    free(mem);
    free(inbuff);
    free(outbuff);
    free(outbuff_gold);

    return 0;
}

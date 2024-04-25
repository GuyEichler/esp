// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "../inc/espacc_config.h"
#include "../inc/espacc.h"
#include "../inc/espacc_tb_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

// #include "A_array.h"
// #include "H_array.h"
// #include "W_array.h"
// #include "Q_array.h"
// #include "initial_state_array.h"
// #include "measurements_array.h"
// #include "prediction_array.h"
// #include "real_array.h"
// #include "P_array.h"

// #include "A_array_soma.h"
// #include "H_array_soma.h"
// #include "W_array_soma.h"
// #include "Q_array_soma.h"
// #include "initial_state_array_soma.h"
// #include "measurements_array_soma.h"
// #include "prediction_array_soma.h"
// #include "real_array_soma.h"
// #include "P_array_soma.h"

#include "A_array_hc.h"
#include "H_array_hc.h"
#include "W_array_hc.h"
#include "Q_array_hc.h"
#include "initial_state_array_hc.h"
#include "measurements_array_hc.h"
#include "prediction_array_hc.h"
#include "real_array_hc.h"
#include "P_array_hc.h"

#include <eigen3/Eigen/Dense>
using namespace Eigen;

// #define STATES 6
// #define NEURONS 164
// #define TIME_STAMPS 100
// #define CHUNKS 100
// #define BATCHES TIME_STAMPS / CHUNKS

// #define STATES 6
// #define NEURONS 52
// #define TIME_STAMPS 100
// #define CHUNKS 100
// #define BATCHES TIME_STAMPS / CHUNKS

#define STATES 6
#define NEURONS 46
#define TIME_STAMPS 100
#define CHUNKS 100
#define BATCHES TIME_STAMPS / CHUNKS

int main(int argc, char **argv) {

    printf("****start*****\n");

    /* <<--params-->> */
    const unsigned inv_reset = 6;
    const unsigned inv_num = 6;
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


    in_words_adj = round_up(x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim, VALUES_PER_WORD);
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
                inbuff[i * in_words_adj + j] = (word_t) A[j - (x_dim + x_dim * x_dim)];
            }

            //Q
            for(; j < x_dim + x_dim * x_dim * 3; j++)
            {
                inbuff[i * in_words_adj + j] = (word_t) W[j - (x_dim + x_dim * x_dim * 2)];
            }

            //R
            for(; j < x_dim + x_dim * x_dim * 3 + z_dim * z_dim; j++)
            {
                inbuff[i * in_words_adj + j] = (word_t) Q[j - (x_dim + x_dim * x_dim * 3)];
            }

            //H
            for(; j < x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim; j++)
            {
                inbuff[i * in_words_adj + j] = (word_t) H[j - (x_dim + x_dim * x_dim * 3 + z_dim * z_dim)];
                //printf("Value of H = %f \n", measurements[NEURONS * (i+1) + j]);
            }
        }

        unsigned base_index = (x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim);
        //Z
        if (i == 0)
            for(; j < base_index + z_dim * chunks; j++)
            {
                inbuff[j] = (word_t) measurements[NEURONS * (i+1) + j - base_index];

                // printf("Value of Z = %f index %d inbuff %d\n", measurements[NEURONS * (i+1) + j - base_index], NEURONS * (i+1) + j - base_index, j);
            }
        else
            for(; j < z_dim * chunks; j++)
            {
                inbuff[in_words_adj + i * in_words_adj_z + j] = (word_t) measurements[NEURONS * (i*chunks+1) + j];

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


    for(int i = 0; i < TIME_STAMPS; i++)
    {
        Matrix<float, NEURONS, NEURONS> Mat_S = Mat_H * (Mat_F * Mat_P * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() + Mat_R;
        Matrix<float, STATES, NEURONS> Mat_K = (Mat_F * Mat_P * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() * Mat_S.inverse();
        Mat_P = (Mat_I -  Mat_K * Mat_H) * (Mat_F * Mat_P * Mat_F.transpose() + Mat_Q);
        // if(i == 0)
        //     std::cout << "Matrix P = " << Mat_P << std::endl;

        float P_flat[STATES*STATES];
        Map<Matrix<float, STATES, STATES, RowMajor> >(&P_flat[0], Mat_P.rows(), Mat_P.cols()) = Mat_P;

        for(unsigned j = 0; j < x_dim + x_dim * x_dim; j++)
        {
            if(j < x_dim)
                outbuff_gold[i * out_words_adj/chunks + j] = (word_t) prediction[STATES * (i+1) + j];
            else
            {
                outbuff_gold[i * out_words_adj/chunks + j] = (word_t) P_flat[j - x_dim];
            }
        }

        // Matrix<float, NEURONS, NEURONS> Mat_S_iter;

        // //Norm check

        // //I - A*X_0
        // Matrix<float, NEURONS, NEURONS> res = Mat_S * Mat_S.inverse();
        // Matrix<float, NEURONS, NEURONS> I_res = Matrix<float, NEURONS, NEURONS>::Identity();
        // Matrix<float, NEURONS, NEURONS> final_res = I_res - res;
        // float norm = final_res.norm();
        // std::cout << "Iteration: " << i << " Norm for calculation = " << norm << std::endl;



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

        // int inv_ok = inverse_tb<Z_MAX, float>(S_array_gauss, S_inv_array_gauss, z_dim);

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
        // std::cout << "Iteration: " << i << " Norm for gauss = " << norm_gauss << std::endl;
        // std::cout << "Iteration: " << i << " max diff from gauss = " << diff_g << std::endl;




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
        //     Matrix<float, NEURONS, NEURONS> diff_mat_S = Mat_S_iter - Mat_S;
        //     float diff = std::max(diff_mat_S.maxCoeff(), (-1)*diff_mat_S.minCoeff());
        //     std::cout << "Iteration: " << i << " Norm for approximation = " << norm_iter << std::endl;
        //     std::cout << "Iteration: " << i << " max diff = " << diff << std::endl;


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
        //             iterative_inverse_tb<Z_MAX>(Mat_S_array, S_inv_array, S_inv_final_array, z_dim);
        //             //std::cout << "Iterative 0" << std::endl;
        //         }
        //         else if(l % 2 == 1){
        //             iterative_inverse_tb<Z_MAX>(Mat_S_array, S_inv_final_array, S_inv_final2_array, z_dim);
        //             //std::cout << "Iterative 1" << std::endl;
        //         }
        //         else{// if(i % 2 == 0){
        //             iterative_inverse_tb<Z_MAX>(Mat_S_array, S_inv_final2_array, S_inv_final_array, z_dim);
        //             //std::cout << "Iterative 2" << std::endl;
        //         }
        //     }

        //     if(inv_num % 2 == 0)
        //     {//final2
        //         Map<Matrix<float, NEURONS, NEURONS, RowMajor> > iterative_out(&S_inv_final2_array[0][0], NEURONS, NEURONS);
        //         prev_S_inv_iter = iterative_out;
        //         //std::cout << "Set final 2" << std::endl;
        //     }
        //     else
        //     {//final
        //         Map<Matrix<float, NEURONS, NEURONS, RowMajor> > iterative_out(&S_inv_final_array[0][0], NEURONS, NEURONS);
        //         prev_S_inv_iter = iterative_out;
        //         //std::cout << "Set final 1" << std::endl;
        //     }

        //     Matrix<float, STATES, NEURONS> Mat_K_iter = (Mat_F * Mat_P_iter * Mat_F.transpose() + Mat_Q) * Mat_H.transpose() * prev_S_inv_iter;//Mat_S.inverse();
        //     Mat_P_iter = (Mat_I -  Mat_K_iter * Mat_H) * (Mat_F * Mat_P_iter * Mat_F.transpose() + Mat_Q);
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

    }

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

            // if(j%(x_dim + x_dim*x_dim) < x_dim)
            //     std::cout << "OKAY: X Accelerator value: " << acc_val << " Golden value: " << gold_val << " index: " << i * out_words_adj + j << " Iter " << i << " diff: " << diff << std::endl;
            // else
            //     std::cout << "P Accelerator value: " << acc_val << " Golden value: " << gold_val << " index: " << i * out_words_adj + j << " Iter " << i << " diff: " << diff << std::endl;

	    if (outbuff[i * out_words_adj + j] != outbuff_gold[i * out_words_adj + j]){
                if(diff/gold_val > 0.1 || diff/acc_val > 0.1 || diff/gold_val < -0.1 || diff/acc_val < -0.1){
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

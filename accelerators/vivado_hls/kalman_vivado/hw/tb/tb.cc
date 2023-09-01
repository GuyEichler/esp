// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "../inc/espacc_config.h"
#include "../inc/espacc.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "A_array.h"
#include "H_array.h"
#include "W_array.h"
#include "Q_array.h"
#include "initial_state_array.h"
#include "measurements_array.h"
#include "prediction_array.h"
#include "real_array.h"

#include <eigen3/Eigen/Dense>
using namespace Eigen;

#define STATES 6
#define NEURONS 164
#define TIME_STAMPS 10

int main(int argc, char **argv) {

    printf("****start*****\n");

    /* <<--params-->> */
    const unsigned iter = TIME_STAMPS;
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


    in_words_adj = round_up(z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim, VALUES_PER_WORD);
    in_words_adj_z = round_up(z_dim, VALUES_PER_WORD);
    out_words_adj = round_up(x_dim + x_dim * x_dim, VALUES_PER_WORD);
    in_size = in_words_adj + in_words_adj_z * (iter - 1);
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

        //Z
        for(; j < z_dim; j++)
        {
            if(i == 0)
                inbuff[i * in_words_adj + j] = (word_t) measurements[NEURONS * (i+1) + j];
            else
                inbuff[in_words_adj + (i-1) * in_words_adj_z + j] = (word_t) measurements[NEURONS * (i+1) + j];
            // if(i == 3)
            //     printf("Value of Z = %f index %d\n", measurements[NEURONS * (i+1) + j], in_words_adj + (i-1) * in_words_adj_z + j);
        }

        if(i == 0) //only for first iteration
        {
            //X
            for(; j < z_dim + x_dim; j++)
            {
                inbuff[i * in_words_adj + j] = (word_t) initial[j - z_dim];
            }

            //P
            for(; j < z_dim + x_dim + x_dim * x_dim; j++)
            {
                inbuff[i * in_words_adj + j] = (word_t) 0.0;
            }

            //F
            for(; j < z_dim + x_dim + x_dim * x_dim * 2; j++)
            {
                inbuff[i * in_words_adj + j] = (word_t) A[j - (z_dim + x_dim + x_dim * x_dim)];
            }

            //Q
            for(; j < z_dim + x_dim + x_dim * x_dim * 3; j++)
            {
                inbuff[i * in_words_adj + j] = (word_t) W[j - (z_dim + x_dim + x_dim * x_dim * 2)];
            }

            //R
            for(; j < z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim; j++)
            {
                inbuff[i * in_words_adj + j] = (word_t) Q[j - (z_dim + x_dim + x_dim * x_dim * 3)];
            }

            //H
            for(; j < z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim; j++)
            {
                inbuff[i * in_words_adj + j] = (word_t) H[j - (z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim)];
                //printf("Value of H = %f \n", measurements[NEURONS * (i+1) + j]);
            }
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


    for(int i = 0; i < iter; i++)
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
                outbuff_gold[i * out_words_adj + j] = (word_t) prediction[STATES * (i+1) + j];
            else
            {
                outbuff_gold[i * out_words_adj + j] = (word_t) P_flat[j - x_dim];
            }
        }

    }

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
    for(unsigned i = 0; i < iter; i++)
        for(unsigned j = 0; j < x_dim + x_dim * x_dim; j++)
        {
            word_t gold_val = outbuff_gold[i * out_words_adj + j];
            word_t acc_val = outbuff[i * out_words_adj + j];
            //printf("Accelerator value: %f Golden value: %f index: %d\n", acc_val, gold_val, i * out_words_adj + j);

            word_t diff = std::abs(gold_val - acc_val);
            MSE += std::pow(diff, 2.0);

            if(j < x_dim)
                std::cout << "X Accelerator value: " << acc_val << " Golden value: " << gold_val << " index: " << i * out_words_adj + j << " Iter " << i << " diff: " << diff << std::endl;
            // else
            //     std::cout << "P Accelerator value: " << acc_val << " Golden value: " << gold_val << " index: " << i * out_words_adj + j << " Iter " << i << " diff: " << diff << std::endl;

	    if (outbuff[i * out_words_adj + j] != outbuff_gold[i * out_words_adj + j]){
                if(diff/gold_val > 0.5){
                    //printf("Accelerator value: %f Golden value: %f index: %d\n", acc_val, gold_val, i * out_words_adj + j);
                    if(j < x_dim)
                        std::cout << "X Accelerator value: " << acc_val << " Golden value: " << gold_val << " index: " << i * out_words_adj + j << " Iter " << i << " diff: " << diff << std::endl;
                    else
                        std::cout << "P Accelerator value: " << acc_val << " Golden value: " << gold_val << " index: " << i * out_words_adj + j << " Iter " << i << " diff: " << diff << std::endl;

                    errors++;
                }
            }
        }

    MSE /= ((x_dim + x_dim * x_dim) * iter);
    std::cout << "Output MSE: " << MSE << std::endl;

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

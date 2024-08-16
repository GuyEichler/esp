// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "../inc/espacc_config.h"
#include "../inc/espacc.h"
#include "../inc/espacc_functions.h"
#include "hls_stream.h"
#include "hls_math.h"
#include <cstring>

void load(word_t Z[CHUNK_MAX][Z_MAX],
          word_t X[1][X_MAX],
          word_t P[X_MAX][X_MAX],
          word_t F[X_MAX][X_MAX],
          word_t Q_kal[X_MAX][X_MAX],
          word_t R_kal[Z_MAX][Z_MAX],
          word_t H[Z_MAX][X_MAX],
          // word_t X_pred[1][X_MAX],
          // word_t P_pred[X_MAX][X_MAX],
          // word_t _inbuff[SIZE_IN_CHUNK_DATA],
          dma_word_t *in1,
          /* <<--compute-params-->> */
          const unsigned chunks,
          const unsigned iter,
          const unsigned x_dim,
          const unsigned z_dim,
	  dma_info_t &load_ctrl,
          unsigned batch,
          // unsigned curr_chunk,
          bool &enable)
{
load_data:

    // unsigned length_raw = z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim;
    unsigned total_length =
        round_up(x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim, VALUES_PER_WORD);

    unsigned z_length = round_up(z_dim * chunks, VALUES_PER_WORD);

    //const unsigned index = 0; //length * (batch * 1 + chunk);
    unsigned length_Z = z_dim * chunks;
    unsigned length_X = x_dim;
    unsigned length_P = x_dim * x_dim;
    unsigned length_F = x_dim * x_dim;
    unsigned length_Q = x_dim * x_dim;
    unsigned length_R = z_dim * z_dim;
    unsigned length_H = z_dim * x_dim;

    // //Check if we're in first iteraion or not
    // unsigned batch_0 = batch == 0 ? 0 : 1;
    // unsigned batch_1 = batch > 1 ? (batch - 1) : 0;
    // unsigned batch_total = batch == 0 ? 1 : 0;

    // unsigned final_index = 0;// = total_length * batch_0;
    // //if(batch == 0) final_index = 0;
    // if (batch > 0) final_index = total_length + z_length * (batch - 1);
    //else if (batch > 1) final_index = total_length + z_length * (batch - 1);
    // unsigned index_z;// = z_length * batch_1;
    // if(batch > 1) index_z = z_length * (batch - 1);
    // else index_z = 0;

    // const unsigned index = final_index;
    unsigned index = (batch == 0) ? 0 : (total_length + z_length * (batch));
    // unsigned length = 0;// = total_length * batch_total + z_length * batch_0;
    // if(batch == 0) length = total_length;
    // else length = z_length;
    unsigned length = (batch == 0) ? total_length + z_length : z_length;

    unsigned X_index = index;//Z_index + length_Z;
    unsigned P_index = X_index + length_X;
    unsigned F_index = P_index + length_P;
    unsigned Q_index = F_index + length_F;
    unsigned R_index = Q_index + length_Q;
    unsigned H_index = R_index + length_R;
    unsigned Z_index = H_index + length_H;// + index_z;
    unsigned end_index = Z_index + length_Z;

    unsigned row = 0;
    word_t tmp[2];
    word_t dummy;
    unsigned base_index = (batch == 0) ? 0 : Z_index;

    unsigned dma_length = length / VALUES_PER_WORD;
    unsigned dma_index = index / VALUES_PER_WORD;

    if(enable){

    // if(curr_chunk == 0){

    load_ctrl.index = dma_index;
    load_ctrl.length = dma_length;
    load_ctrl.size = SIZE_WORD_T;

#ifndef __SYNTHESIS__
    printf("START LOADING \n");
    printf("dma_length %d\n", dma_length);
    printf("dma_index %d\n", dma_index);
#endif

    if(batch == 0)
    {

// #ifndef __SYNTHESIS__
//     printf("INIT ARRAYS \n");
// #endif
    load_init0:for(unsigned k = 0; k < Z_MAX; k++)
    load_init1:for(unsigned j = 0; j < Z_MAX; j++)
        {

            if(j == 0)
            {
                if(k < X_MAX)
                {
                    X[0][k] = 0.0;
                    //X_pred[0][k] = 0.0;
                }
            }

            if(j < chunks)
                Z[j][k] = 0.0;

            if(k < X_MAX && j < X_MAX)
            {
                P[k][j] = 0.0;
                //P_pred[k][j] = 0.0;
                F[k][j] = 0.0;
                Q_kal[k][j] = 0.0;
            }

            if (j < X_MAX)
                H[k][j] = 0.0;

            R_kal[k][j] = 0.0;
        }

    }

// #ifndef __SYNTHESIS__
//     printf("DMA LOAD \n");
// #endif

    for (unsigned i = 0; i < dma_length; i++) {
#pragma HLS loop_tripcount max=22197

    load_label0:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
            // _inbuff[i * VALUES_PER_WORD + j] = in1[dma_index + i].word[j];

            tmp[j] = in1[dma_index + i].word[j];
            unsigned col_index = 0;
            unsigned total_index = (i * VALUES_PER_WORD) + j + base_index;

            //Make the arrays 2D with an offset\counter that jumps when needed
//             if (total_index < X_index){//Z
//                 col_index = total_index;
//                 row = 0;
//                 Z[0][col_index] = tmp[j];
// // #ifndef __SYNTHESIS__
// //                 printf("Loaded Z: %f index %d\n", tmp[j], i * VALUES_PER_WORD + j);
// // #endif
//             }
            if (total_index < P_index){//X
                col_index = total_index - X_index;
                row = 0;
                X[0][col_index] = tmp[j];
            }
            else if (total_index < F_index){//P
                col_index = total_index - P_index - row*x_dim;
                P[row][col_index] = tmp[j];
                if(col_index == x_dim - 1)
                    row++;
                if(row == x_dim)
                    row = 0;
            }
            else if (total_index < Q_index){//F
                col_index = total_index - F_index - row*x_dim;
                F[row][col_index] = tmp[j];
                if(col_index == x_dim - 1)
                    row++;
                if(row == x_dim)
                    row = 0;
            }
            else if (total_index < R_index){//Q
                col_index = total_index - Q_index - row*x_dim;
                Q_kal[row][col_index] = tmp[j];
                if(col_index == x_dim - 1)
                    row++;
                if(row == x_dim)
                    row = 0;
            }
            else if (total_index < H_index){//R
                col_index = total_index - R_index - row*z_dim;
                R_kal[row][col_index] = tmp[j];
                if(col_index == z_dim - 1)
                    row++;
                if(row == z_dim)
                    row = 0;
            }
            else if (total_index < Z_index){//H
                col_index = total_index - H_index - row*x_dim;
                H[row][col_index] = tmp[j];
                if(col_index == x_dim - 1)
                    row++;
                if(row == z_dim)
                    row = 0;
            }
            else if (total_index < end_index){//Z
                col_index = total_index - Z_index - row*z_dim;
                Z[row][col_index] = tmp[j];
// #ifndef __SYNTHESIS__
//                 if(batch == 1)
//                 printf("Loaded Z[%d][%d]: %f\n", row, col_index, Z[row][col_index]);
// #endif
                if(col_index == z_dim - 1)
                    row++;
                if(row == chunks)
                    row = 0;
            }
            else{
                dummy = tmp[j];
                // printf("index read by tmp is: %d, dma_length: %d \n", total_index, i+1);
            }
    	}
    }
    }

//     if(batch != 0 || curr_chunk != 0)
//     {
// // #ifndef __SYNTHESIS__
// //     printf("COPY PREDICTIONS \n");
// // #endif
//         //Copy X/P_pred into X/P in case we have more than one iteration
//         row = 0;
//         load_pred:for(unsigned i = 0; i < (x_dim + x_dim * x_dim); i++)
//         {
// #pragma HLS loop_tripcount max=42

//             word_t tmp_pred;
//             if(i < x_dim)
//             {
//                 tmp_pred = X_pred[0][i];
//                 X[0][i] = tmp_pred;

// // #ifndef __SYNTHESIS__
// //                 printf("Copied value X_pred: %f i=%d\n", tmp_pred, i);
// // #endif

//             }
//             else
//             {
//                 unsigned col = i - x_dim - row*x_dim;
//                 tmp_pred = P_pred[row][col];
//                 P[row][col] = tmp_pred;

//                 if(col == x_dim - 1)
//                     row++;

// // #ifndef __SYNTHESIS__
// //                 printf("Copied value P_pred: %f \n", tmp);
// // #endif
//             }
//         }
//     }
//}

}

void store(// word_t X_pred[1][X_MAX],
           // word_t P_pred[X_MAX][X_MAX],
           word_t _outbuff[SIZE_OUT_CHUNK_DATA],
           dma_word_t *out,
          /* <<--compute-params-->> */
           const unsigned chunks,
           const unsigned iter,
           const unsigned x_dim,
           const unsigned z_dim,
	   dma_info_t &store_ctrl,
           unsigned batch,
           // unsigned curr_chunk,
           bool &enable)
{
store_data:

    const unsigned length =
        round_up((x_dim + x_dim * x_dim) * chunks, VALUES_PER_WORD);
    const unsigned store_offset =
        round_up(x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim , VALUES_PER_WORD) +
        round_up(z_dim*chunks, VALUES_PER_WORD) * iter;
    const unsigned out_offset = store_offset;
    const unsigned index = out_offset + length * batch;

    unsigned dma_length = length / VALUES_PER_WORD;
    unsigned dma_index = index / VALUES_PER_WORD;

    if(enable){// && (curr_chunk == chunks-1)){

#ifndef __SYNTHESIS__
    printf("START STORING \n");
    printf("dma_length %d\n", dma_length);
    printf("dma_index %d\n", dma_index);
#endif

    store_ctrl.index = dma_index;
    store_ctrl.length = dma_length;
    store_ctrl.size = SIZE_WORD_T;

    // int row = 0;

    for (unsigned i = 0; i < dma_length; i++) {
#pragma HLS loop_tripcount max=2100

    store_label1:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
	    out[dma_index + i].word[j] = _outbuff[i * VALUES_PER_WORD + j];

// #ifndef __SYNTHESIS__
//             int index = i * VALUES_PER_WORD + j;

//             if(index%(x_dim + x_dim * x_dim) < x_dim)
//             {
// //                 word_t tmp = X_pred[0][index];
// //                 out[dma_index + i].word[j] = tmp;


//                 printf("outbuff[%d] = %f\n", index, _outbuff[index]);

//             }
// #endif
//             else
//             {
//                 int col = index - x_dim - row*x_dim;
//                 word_t tmp = P_pred[row][col];
//                 out[dma_index + i].word[j] = tmp;

// #ifndef __SYNTHESIS__
//                 // printf("outbuff[%d] = %.12f , P_pred[%d][%d] = %.12f \n", i, _outbuff[i], row, col, P_pred[row][col]);
// #endif

//                 if(col == x_dim - 1)
//                     row++;

//             }


// #ifndef __SYNTHESIS__
//             // printf("out[%d].word[%d] = %f\n", i+dma_index, j, _outbuff[i * VALUES_PER_WORD + j]);
// #endif
	}
    }
    }
}


void compute(word_t Z[CHUNK_MAX][Z_MAX],
             word_t X[1][X_MAX],
             word_t P[X_MAX][X_MAX],
             word_t F[X_MAX][X_MAX],
             word_t Q_kal[X_MAX][X_MAX],
             word_t R_kal[Z_MAX][Z_MAX],
             word_t H[Z_MAX][X_MAX],
             //word_t _inbuff[SIZE_IN_CHUNK_DATA],
             /* <<--compute-params-->> */
             const unsigned inv_reset,
             const unsigned inv_num,
             const unsigned chunks,
             const unsigned iter,
             const unsigned x_dim,
             const unsigned z_dim,
             // word_t X_pred[1][X_MAX],
             // word_t P_pred[X_MAX][X_MAX],
             word_t _outbuff[SIZE_OUT_CHUNK_DATA],
             // unsigned curr_chunk,
             unsigned curr_batch,
             bool &enable)
{

    // TODO implement compute functionality
    // const unsigned length =
    //     round_up(z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim, VALUES_PER_WORD) / 1;

    // for (int i = 0; i < length; i++)
    //     _outbuff[i] = _inbuff[i];

    // unsigned XX_size = x_dim*x_dim;
    // unsigned ZZ_size = z_dim*z_dim;
    // unsigned ZX_size = z_dim*x_dim;

    word_t inter1[X_MAX][X_MAX];//Define inter1 with the size of inter3 to be
    //able to share resources
    // word_t inter1[Z_MAX][X_MAX];
    word_t inter2[X_MAX][X_MAX];

    word_t inter3[Z_MAX][X_MAX];//Reuse inter 1
    word_t S[Z_MAX][Z_MAX];

    static word_t S_inv[Z_MAX][Z_MAX];
    static word_t S_inv_final[Z_MAX][Z_MAX];
    static word_t S_inv_final2[Z_MAX][Z_MAX];
    int inv_ok = 0;

    //Compute K = inter2 * H.T * S_inv
    word_t inter4[X_MAX][Z_MAX];//Reuse inter3 by changing computation
    word_t K[X_MAX][Z_MAX];

    //Compute Y = Z - (H * F * X)
    word_t inter5[Z_MAX][X_MAX];//Reuse inter3
    word_t Y[Z_MAX][1];

    //Compute final prediction state X_pred = F * X + K * Y
    word_t inter6[1][X_MAX];//Reuse inter1

    //Compute final prediction state covariance P_pred = (I - K * H) * inter2
    //static word_t P_pred[X_MAX][X_MAX];//Set as PLMs from top
    word_t inter7[X_MAX][X_MAX];// Reuse inter1

    static word_t X_pred_ping[1][X_MAX];
    static word_t X_pred_pong[1][X_MAX];
    static word_t P_pred_ping[X_MAX][X_MAX];
    static word_t P_pred_pong[X_MAX][X_MAX];
    static bool pingpong = true;

    static unsigned inv_counter = 1;

    unsigned in_length = x_dim + x_dim * x_dim;
    if(curr_batch == 0)
    init_pred_1:for(unsigned i = 0; i < X_MAX; i++)
    init_pred_2:for(unsigned j = 0; j < X_MAX; j++)
    {
        word_t tmp_pred_x;
        word_t tmp_pred_p;
        if(i == 0 && j < x_dim){
            tmp_pred_x = X[0][j];
            if(pingpong){
                X_pred_ping[0][j] = tmp_pred_x;
                X_pred_pong[0][j] = 0.0;
            }else{
                X_pred_pong[0][j] = tmp_pred_x;
                X_pred_ping[0][j] = 0.0;
            }
        }
        else if(i == 0){
            X_pred_ping[0][j] = 0.0;
            X_pred_pong[0][j] = 0.0;
        }

        if(i < x_dim && j < x_dim)
        {
            tmp_pred_p = P[i][j];
            if(pingpong){
                P_pred_ping[i][j] = tmp_pred_p;
                P_pred_pong[i][j] = 0.0;
            }else{
                P_pred_pong[i][j] = tmp_pred_p;
                P_pred_ping[i][j] = 0.0;
            }
        }
        else{
            P_pred_pong[i][j] = 0.0;
            P_pred_ping[i][j] = 0.0;
        }
    }

    if(curr_batch == 0)
    {
        inv_counter = 1;
    }

    compute_loop:

    for(int c = 0; c < CHUNK_MAX; c++){

    if(c < chunks){
    // if(enable){

#ifndef __SYNTHESIS__
    printf("inv_counter = %d\n", inv_counter);
    printf("curr_batch = %d\n", curr_batch);
    printf("c = %d\n", c);
    printf("pingpong = %d\n", pingpong);
#endif

// #ifndef __SYNTHESIS__
//     printf("Z = \n");
//     hls::print_matrix<CHUNK_MAX, Z_MAX, word_t, hls::NoTranspose>((word_t(*)[Z_MAX])Z, "   ");
// #endif

    //Compute inter1 = F x P
    if(pingpong)
    mul_int11:hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             X_MAX, X_MAX, X_MAX, X_MAX,
                             X_MAX, X_MAX,
                             TRAITS_A,
                             word_t, word_t> (F, P_pred_ping, inter1);
    else
    mul_int12:hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             X_MAX, X_MAX, X_MAX, X_MAX,
                             X_MAX, X_MAX,
                             TRAITS_A,
                             word_t, word_t> (F, P_pred_pong, inter1);

#ifndef __SYNTHESIS__
        // printf("inter1 = \n");
        // hls::print_matrix<X_MAX, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])inter1, "   ");
#endif

    //Compute inter2 = inter1 x F.T
    hls::matrix_multiply_top<hls::NoTranspose, hls::Transpose,
                             X_MAX, X_MAX, X_MAX, X_MAX,
                             X_MAX, X_MAX,
                             TRAITS_B,
                             word_t, word_t> (inter1, F, inter2);

#ifndef __SYNTHESIS__
    // printf("inter2 = \n");
    // hls::print_matrix<X_MAX, X_MAX, word_t,
    // hls::NoTranspose>((word_t(*)[X_MAX])inter2, "   ");
    // printf("F = \n");
    // hls::print_matrix<X_MAX, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])F, "   ");
#endif

    //Compute inter2 = inter2 + Q
    LOOP_INT2_1:for(int i = 0; i < X_MAX; i++)
    LOOP_INT2_2:for(int j = 0; j < X_MAX; j++)
        {
            if(i < x_dim && j < x_dim)
            {
                word_t tmp = inter2[i][j];
                inter2[i][j] = tmp + Q_kal[i][j];
            }
            else
                inter2[i][j] = 0.0;
        }

    // compute_int2(inter2, Q_kal, x_dim);

#ifndef __SYNTHESIS__
        // printf("inter2 = \n");
        // hls::print_matrix<X_MAX, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])inter2, "   ");
#endif

    //Compute S = H * inter2 * H.T + R
    // word_t inter3[Z_MAX][X_MAX];//Reuse inter 1
    // word_t S[Z_MAX][Z_MAX];

    //Compute inter3 = H x inter2
    hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             Z_MAX, X_MAX, X_MAX, X_MAX,
                             Z_MAX, X_MAX,
                             TRAITS_C,
                             word_t, word_t> (H, inter2, inter3);

#ifndef __SYNTHESIS__
        // printf("inter3 = \n");
        // hls::print_matrix<Z_MAX, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])inter3, "   ");
#endif

    //Compute S = inter3 * H.T
    hls::matrix_multiply_top<hls::NoTranspose, hls::Transpose,
                             Z_MAX, X_MAX, Z_MAX, X_MAX,
                             Z_MAX, Z_MAX,
                             TRAITS_D,
                             word_t, word_t> (inter3, H, S);

    //Compute S = S + R
    LOOP_SR_1:for(int i = 0; i < Z_MAX; i++)
    LOOP_SR_2:for(int j = 0; j < Z_MAX; j++)
        {
            if(i < z_dim && j < z_dim)
            {
                word_t tmp = S[i][j];
                S[i][j] = tmp + R_kal[i][j];
            }
            // else if(i == j)//(i >= z_dim )
            // {
            //     S[i][j] = 1.0;
            // }
            else
                S[i][j] = 0.0;
        }
    // compute_sr(S, R_kal, z_dim);

    //Compute S^-1 = S_inv - QR inverse
    // word_t S_inv[Z_MAX][Z_MAX];
    // int inv_ok = 0;

    //Initialize S_inv as Identity
    // if((curr_batch == 0 && c == 0) || inv_num == 0)
    LOOP_S_inv_1:for(int i = 0; i < Z_MAX; i++)
    LOOP_S_inv_2:for(int j = 0; j < Z_MAX; j++)
        {
            if(i == j && i < z_dim)
            {
                if((curr_batch == 0 && c == 0) || inv_num == 0 || inv_counter == inv_reset)
                {
                    S_inv[i][j] = 1.0;
                    S_inv_final[i][j] = 1.0;
                    S_inv_final2[i][j] = 1.0;
                    // after_reset = true;
                }
                //else if(after_reset == false)
                else if ((inv_counter > 1 || inv_reset == 0) && (curr_batch * chunks + c != 1))
                //else if ((curr_batch != 0 || c > 1) && (inv_counter > 1 || inv_reset == 0))
                {//Comment this block to change policy
                    //printf("inv_counter is: %d\n", inv_counter);
                    word_t tmp = S_inv_final[i][j];
                    word_t tmp2 = S_inv_final2[i][j];
                    if(inv_num % 2 == 0)// && c > 1)
                        S_inv[i][j] = tmp2;
                    else// if (c > 1)
                        S_inv[i][j] = tmp;
                }
            }
            else
            {
                if((curr_batch == 0 && c == 0) || inv_num == 0 || inv_counter == inv_reset)
                {
                    S_inv[i][j] = 0.0;
                    S_inv_final[i][j] = 0.0;
                    S_inv_final2[i][j] = 0.0;
                }
                //else if(after_reset == false)
                else if ((inv_counter > 1 || inv_reset == 0) && (curr_batch * chunks + c != 1))
                //else if ((curr_batch != 0 || c > 1) && (inv_counter > 1 || inv_reset == 0))
                {//Comment this block to change policy
                    word_t tmp = S_inv_final[i][j];
                    word_t tmp2 = S_inv_final2[i][j];
                    if(inv_num % 2 == 0)// && c > 1)
                        S_inv[i][j] = tmp2;
                    else// if (c > 1)
                        S_inv[i][j] = tmp;
                }

                // if(i == Z_MAX && j == Z_MAX && after_reset == true)
                //     after_reset = false;
            }

        }
    // init_s_inv(S_inv, z_dim);

// #ifndef __SYNTHESIS__
//     printf("S_inv = \n");
//     hls::print_matrix<Z_MAX, Z_MAX, word_t, hls::NoTranspose>((word_t(*)[Z_MAX])S_inv, "   ");
// #endif

    if((curr_batch == 0 && c == 0) || inv_num == 0 || inv_counter == inv_reset)
    {
        inv_ok = inverse<Z_MAX, word_t>(S, S_inv, z_dim);
        // hls::cholesky_inverse_top<Z_MAX,MY_CONFIG_C,word_t,word_t>(S,S_inv,inv_ok);
        //Use QR decomposition to compute inverse
        //hls::qr_inverse_top<Z_MAX,MY_CONFIG,word_t,word_t>(S,S_inv,inv_ok);
        //hls::qr_inverse<Z_MAX,word_t,word_t>(S,S_inv,inv_ok);
    }else{
        // inverse<Z_MAX, word_t>(S, S_inv, z_dim, inv_ok);
        // inv_ok = inverse_partial<Z_MAX, word_t>(S, S_inv, z_dim);
        //many_iterative_inverse<Z_MAX,1>(S, S_inv, z_dim);
        for(unsigned i = 0; i < inv_num; i++){
        #pragma HLS loop_tripcount max=3
            if(i == 0){
                iterative_inverse<Z_MAX>(S, S_inv, S_inv_final, z_dim);
            }
            else if(i % 2 == 1){
                iterative_inverse<Z_MAX>(S, S_inv_final, S_inv_final2, z_dim);
            }
            else{// if(i % 2 == 0){
                iterative_inverse<Z_MAX>(S, S_inv_final2, S_inv_final, z_dim);
            }
        }

    // iterative_inverse3<Z_MAX>(S, S_inv, S_inv_final, S_inv_final2, z_dim);
    // iterative_inverse<Z_MAX>(S, S_inv_final2, S_inv_final, z_dim);
    // iterative_inverse<Z_MAX>(S, S_inv_final, S_inv_final2, z_dim);
    // iterative_inverse<Z_MAX>(S, S_inv_final2, S_inv_final, z_dim);
    }


#ifndef __SYNTHESIS__
    //printf("inv_ok = %d\n", inv_ok);
    // printf("S_inv = \n");
    // hls::print_matrix<Z_MAX, Z_MAX, word_t, hls::NoTranspose>((word_t(*)[Z_MAX])S_inv, "   ");
#endif

    // //Compute K = inter2 * H.T * S_inv
    // word_t inter4[X_MAX][Z_MAX];//Reuse inter3 by changing computation
    // word_t K[X_MAX][Z_MAX];

    //Compute inter4 = inter2 * H.T
    hls::matrix_multiply_top<hls::NoTranspose, hls::Transpose,
                             X_MAX, X_MAX, Z_MAX, X_MAX,
                             X_MAX, Z_MAX,
                             TRAITS_J,
                             word_t, word_t> (inter2, H, inter4);
    // //Compute inter4.T = inter3 = H * inter2.T
    // hls::matrix_multiply_top<hls::NoTranspose, hls::Transpose,
    //                          Z_MAX, X_MAX, X_MAX, X_MAX,
    //                          Z_MAX, X_MAX,
    //                          TRAITS_J,
    //                          word_t, word_t> (H, inter2, inter3);

#ifndef __SYNTHESIS__
    // printf("inter4 = \n");
    // hls::print_matrix<X_MAX, Z_MAX, word_t, hls::NoTranspose>((word_t(*)[Z_MAX])inter4, "   ");
#endif

    //Compute K = inter4 * S_inv
    if((curr_batch == 0 && c == 0) || inv_num == 0 || inv_counter == inv_reset)
    hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             X_MAX, Z_MAX, Z_MAX, Z_MAX,
                             X_MAX, Z_MAX,
                             TRAITS_E,
                             word_t, word_t> (inter4, S_inv, K);
                             // word_t, word_t> (inter4, S_inv_final2, K);
    else{
        if(inv_num % 2 == 0){
            hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             X_MAX, Z_MAX, Z_MAX, Z_MAX,
                             X_MAX, Z_MAX,
                             TRAITS_E,
                             // word_t, word_t> (inter4, S_inv, K);
                             word_t, word_t> (inter4, S_inv_final2, K);
        }
        else{
            hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             X_MAX, Z_MAX, Z_MAX, Z_MAX,
                             X_MAX, Z_MAX,
                             TRAITS_E,
                             // word_t, word_t> (inter4, S_inv, K);
                             word_t, word_t> (inter4, S_inv_final, K);
        }
    }


    // mat_mul<X_MAX,Z_MAX>(inter4, S_inv, K, x_dim, z_dim);

    // //Compute K = inter4 * S_inv = inter3.T * S_inv
    // hls::matrix_multiply_top<hls::Transpose, hls::NoTranspose,
    //                          Z_MAX, X_MAX, Z_MAX, Z_MAX,
    //                          X_MAX, Z_MAX,
    //                          TRAITS_E,
    //                          word_t, word_t> (inter3, S_inv, K);

#ifndef __SYNTHESIS__
    // printf("K = \n");
    // hls::print_matrix<X_MAX, Z_MAX, word_t, hls::NoTranspose>((word_t(*)[Z_MAX])K, "   ");
#endif

    // //Compute Y = Z - (H * F * X)
    // word_t inter5[Z_MAX][X_MAX];//Reuse inter3
    // word_t Y[Z_MAX][1];

    //Compute inter5 = H * F
    hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             Z_MAX, X_MAX, X_MAX, X_MAX,
                             Z_MAX, X_MAX,
                             TRAITS_C,
                             word_t, word_t> (H, F, inter5);

    //Compute Y = inter5 * X
    if(pingpong)
    hls::matrix_multiply_top<hls::NoTranspose, hls::Transpose,
                             Z_MAX, X_MAX, 1, X_MAX,
                             Z_MAX, 1,
                             TRAITS_F,
                             word_t, word_t> (inter5, X_pred_ping, Y);
    else
    hls::matrix_multiply_top<hls::NoTranspose, hls::Transpose,
                             Z_MAX, X_MAX, 1, X_MAX,
                             Z_MAX, 1,
                             TRAITS_F,
                             word_t, word_t> (inter5, X_pred_pong, Y);

    //Compute Y = Z - Y
    LOOP_Y:for(int i = 0; i < Z_MAX; i++)
    {
        if(i < z_dim)
        {
            word_t tmp = Y[i][0];
            Y[i][0] = Z[c][i] - tmp;
// #ifndef __SYNTHESIS__
//                 printf("Computed with Z[%d][%d]: %f\n", curr_chunk, i, Z[curr_chunk][i]);
//                 printf("Y[%d][%d]: %f\n", i, 0, Y[i][0]);
// #endif
        }
        else
            Y[i][0] = 0.0;
    }
    // compute_y(Y, Z, z_dim, curr_chunk);

    // //Compute final prediction state X_pred = F * X + K * Y
    // //static word_t X_pred[1][X_MAX];//Set as PLMs from top
    // word_t inter6[1][X_MAX];//Reuse inter1

    //Compute inter6 = K * Y
    hls::matrix_multiply_top<hls::Transpose, hls::Transpose,
                             Z_MAX, 1, X_MAX, Z_MAX,
                             1, X_MAX,
                             TRAITS_G,
                             word_t, word_t> (Y, K, inter6);

    //Compute X_pred = F * X (Compute X.T * F.T)
    if(pingpong)
    hls::matrix_multiply_top<hls::NoTranspose, hls::Transpose,
                             1, X_MAX, X_MAX, X_MAX,
                             1, X_MAX,
                             TRAITS_H,
                             word_t, word_t> (X_pred_ping, F, X_pred_pong);
    else
    hls::matrix_multiply_top<hls::NoTranspose, hls::Transpose,
                             1, X_MAX, X_MAX, X_MAX,
                             1, X_MAX,
                             TRAITS_H,
                             word_t, word_t> (X_pred_pong, F, X_pred_ping);

    //Compute X_pred = X_pred + inter6
    LOOP_X_PRED:for(int i = 0; i < X_MAX; i++)
    {
        if(i < x_dim)
        {
            // word_t tmp = X_pred[0][i];
            word_t tmp2 = inter6[0][i];
            if(pingpong)
                X_pred_pong[0][i] = tmp2 + X_pred_pong[0][i];
            else
                X_pred_ping[0][i] = tmp2 + X_pred_ping[0][i];
// #ifndef __SYNTHESIS__
//                 printf("Compute value X_pred: %f \n", X_pred[0][i]);
// #endif
        }
        else
            if(pingpong)
                X_pred_pong[0][i] = 0.0;
            else
                X_pred_ping[0][i] = 0.0;
    }
    // compute_x_pred(X_pred, inter6, x_dim);

#ifndef __SYNTHESIS__
    // printf("X_pred = \n");
    // hls::print_matrix<1, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])X_pred, "   ");
#endif

    // //Compute final prediction state covariance P_pred = (I - K * H) * inter2
    // //static word_t P_pred[X_MAX][X_MAX];//Set as PLMs from top
    // word_t inter7[X_MAX][X_MAX];// Reuse inter1

    //Compute inter7 = K * H
    hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             X_MAX, Z_MAX, Z_MAX, X_MAX,
                             X_MAX, X_MAX,
                             TRAITS_I,
                             word_t, word_t> (K, H, inter7);

#ifndef __SYNTHESIS__
    // printf("inter7 = \n");
    // hls::print_matrix<X_MAX, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])inter7, "   ");
#endif

    //Compute inter7 = I - inter7
    LOOP_INT7_1:for(int i = 0; i < X_MAX; i++)
    LOOP_INT7_2:for(int j = 0; j < X_MAX; j++)
        {
            if(i < x_dim && j < x_dim)
            {
                word_t tmp = inter7[i][j];
                if(i == j)
                    inter7[i][j] = 1 - tmp;
                else
                    inter7[i][j] = 0 - tmp;
            }
            else
                inter7[i][j] = 0;
        }
    // compute_int7(inter7, x_dim);

    //Compute P_pred = inter7 * inter2
    if(pingpong)
    hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             X_MAX, X_MAX, X_MAX, X_MAX,
                             X_MAX, X_MAX,
                             TRAITS_A,
                             word_t, word_t> (inter7, inter2, P_pred_pong);
    else
    hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             X_MAX, X_MAX, X_MAX, X_MAX,
                             X_MAX, X_MAX,
                             TRAITS_A,
                             word_t, word_t> (inter7, inter2, P_pred_ping);

#ifndef __SYNTHESIS__
    // printf("P_pred = \n");
    // hls::print_matrix<X_MAX, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])P_pred, "   ");
#endif

    //Load output into outbuff
    int out_length = x_dim + x_dim * x_dim;
    int row = 0;

    LOOP_OUT:for(int i = 0; i < X_MAX + X_MAX * X_MAX; i++)
    {
//#pragma HLS loop_tripcount max=42
        if(i < x_dim)
        {
            if(pingpong)
                _outbuff[i + c * out_length] = X_pred_pong[0][i];
            else
                _outbuff[i + c * out_length] = X_pred_ping[0][i];

// #ifndef __SYNTHESIS__
//             printf("compute outbuff[%d] = %f\n", i, _outbuff[i + curr_chunk * out_length]);
// #endif
        }
        else if(i < out_length)
        {
            int col = i - x_dim - row*x_dim;
            if(pingpong)
                _outbuff[i + c * out_length] = P_pred_pong[row][col];
            else
                _outbuff[i + c * out_length] = P_pred_ping[row][col];

#ifndef __SYNTHESIS__
            // printf("outbuff[%d] = %.12f , P_pred[%d][%d] = %.12f \n", i, _outbuff[i], row, col, P_pred[row][col]);
#endif

            if(col == x_dim - 1)
                row++;

        }
    }

#ifndef __SYNTHESIS__
    // printf("X_pred = \n");
    // if(pingpong)
    //     hls::print_matrix<1, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])X_pred_pong, "   ");
    // else
    //     hls::print_matrix<1, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])X_pred_ping, "   ");
#endif

    // compute_out(_outbuff, X_pred, P_pred, x_dim, curr_chunk);
    pingpong = !pingpong;
    if(inv_counter < inv_reset)
        inv_counter++;
    else
        inv_counter = 1;

    }
    }

}


void top(dma_word_t *out, dma_word_t *in1,
         /* <<--params-->> */
	 const unsigned conf_info_inv_reset,
	 const unsigned conf_info_inv_num,
	 const unsigned conf_info_chunks,
	 const unsigned conf_info_iter,
	 const unsigned conf_info_x_dim,
	 const unsigned conf_info_z_dim,
	 dma_info_t &load_ctrl, dma_info_t &store_ctrl)
{

    /* <<--local-params-->> */
	 const unsigned inv_reset = conf_info_inv_reset;
	 const unsigned inv_num = conf_info_inv_num;
	 const unsigned chunks = conf_info_chunks;
	 const unsigned iter = conf_info_iter;
	 const unsigned x_dim = conf_info_x_dim;
	 const unsigned z_dim = conf_info_z_dim;

         //static word_t _inbuff[SIZE_IN_CHUNK_DATA];

         //static word_t Z[1][Z_MAX];
         static word_t Z[CHUNK_MAX][Z_MAX];
         static word_t X[1][X_MAX];
         static word_t P[X_MAX][X_MAX];
         static word_t F[X_MAX][X_MAX];
         static word_t Q_kal[X_MAX][X_MAX];
         static word_t R_kal[Z_MAX][Z_MAX];
         static word_t H[Z_MAX][X_MAX];

         static word_t _outbuff[SIZE_OUT_CHUNK_DATA];

         // static word_t X_pred[1][X_MAX];
         // static word_t P_pred[X_MAX][X_MAX];

         bool enable = false;

    // Batching
batching:
    for (unsigned b = 0; b < iter; b++)
    {
#pragma HLS loop_tripcount max=1
        // Chunking
    go:
        // for (int c = 0; c < chunks; c++)
        // {
//#pragma HLS loop_tripcount max=1
        unsigned curr_b = b;
        // unsigned curr_c = c;

        if(curr_b == 0)
            enable = true;

        load(Z,
             X,
             P,
             F,
             Q_kal,
             R_kal,
             H,
             // X_pred,
             // P_pred,
             //_inbuff,
             in1,
             /* <<--args-->> */
             chunks,
             iter,
             x_dim,
             z_dim,
             load_ctrl,
             curr_b,
             // curr_c,
             enable);

        // if(b == 0)
        compute(Z,
                X,
                P,
                F,
                Q_kal,
                R_kal,
                H,
                //_inbuff,
                /* <<--args-->> */
                inv_reset,
                inv_num,
                chunks,
                iter,
                x_dim,
                z_dim,
                // X_pred,
                // P_pred,
                _outbuff,
                // curr_c,
                curr_b,
                enable);
        // else
        //     compute(Z,
        //             X_pred,
        //             P_pred,
        //             F,
        //             Q,
        //             R,
        //             H,
        //             //_inbuff,
        //             /* <<--args-->> */
        //             iter,
        //             x_dim,
        //             z_dim,
        //             // X_pred,
        //             // P_pred,
        //             _outbuff);

        store(// X_pred,
              // P_pred,
              _outbuff,
              out,
              /* <<--args-->> */
              chunks,
              iter,
              x_dim,
              z_dim,
              store_ctrl,
              curr_b,
              // curr_c,
              enable);

        // }
    }
}


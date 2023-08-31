// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "../inc/espacc_config.h"
#include "../inc/espacc.h"
#include "../inc/espacc_functions.h"
#include "hls_stream.h"
#include "hls_math.h"
#include <cstring>

void load(word_t Z[1][Z_MAX],
          word_t X[1][X_MAX],
          word_t P[X_MAX][X_MAX],
          word_t F[X_MAX][X_MAX],
          word_t Q_kal[X_MAX][X_MAX],
          word_t R_kal[Z_MAX][Z_MAX],
          word_t H[Z_MAX][X_MAX],
          word_t X_pred[1][X_MAX],
          word_t P_pred[X_MAX][X_MAX],
          // word_t _inbuff[SIZE_IN_CHUNK_DATA],
          dma_word_t *in1,
          /* <<--compute-params-->> */
          const unsigned iter,
          const unsigned x_dim,
          const unsigned z_dim,
	  dma_info_t &load_ctrl, int batch)
{
load_data:

    unsigned total_length =
        round_up(z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim, VALUES_PER_WORD);

    unsigned z_length = round_up(z_dim, VALUES_PER_WORD);

    //const unsigned index = 0; //length * (batch * 1 + chunk);
    const unsigned length_Z = z_dim;
    const unsigned length_X = x_dim;
    const unsigned length_P = x_dim * x_dim;
    const unsigned length_F = x_dim * x_dim;
    const unsigned length_Q = x_dim * x_dim;
    const unsigned length_R = z_dim * z_dim;
    const unsigned length_H = z_dim * x_dim;

    //Check if we're in first iteraion or not
    unsigned batch_0 = batch == 0 ? 0 : 1;
    unsigned batch_1 = batch > 1 ? (batch - 1) : 0;
    unsigned batch_total = batch == 0 ? 1 : 0;

    const unsigned index = total_length * batch_0;
    const unsigned index_z = z_length * batch_1;
    const unsigned length = total_length * batch_total + z_length * batch_0;

    unsigned Z_index = index + index_z;
    unsigned X_index = Z_index + length_Z;
    unsigned P_index = X_index + length_X;
    unsigned F_index = P_index + length_P;
    unsigned Q_index = F_index + length_F;
    unsigned R_index = Q_index + length_Q;
    unsigned H_index = R_index + length_R;
    unsigned end_index = H_index + length_H;

    unsigned row = 0;
    word_t tmp[2];
    word_t dummy;

    unsigned dma_length = length / VALUES_PER_WORD;
    unsigned dma_index = (index + index_z) / VALUES_PER_WORD;

    load_ctrl.index = dma_index;
    load_ctrl.length = dma_length;
    load_ctrl.size = SIZE_WORD_T;

#ifndef __SYNTHESIS__
    printf("START LOADING \n");
    printf("dma_length %d\n", dma_length);
    printf("dma_index %d\n", dma_index);
#endif

    for (unsigned i = 0; i < dma_length; i++) {
#pragma HLS loop_tripcount max=14079

    load_label0:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
            // _inbuff[i * VALUES_PER_WORD + j] = in1[dma_index + i].word[j];

            tmp[j] = in1[dma_index + i].word[j];
            unsigned col_index = 0;
            unsigned total_index = (i * VALUES_PER_WORD) + j;

            //Make the arrays 2D with an offset\counter that jumps when needed
            if (total_index < X_index){//Z
                col_index = total_index;
                row = 0;
                Z[0][col_index] = tmp[j];
// #ifndef __SYNTHESIS__
//                 printf("Loaded Z: %f index %d\n", tmp[j], i * VALUES_PER_WORD + j);
// #endif
            }
            else if (total_index < P_index){//X
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
            else if (total_index < end_index){//H
                col_index = total_index - H_index - row*x_dim;
                H[row][col_index] = tmp[j];
                if(col_index == x_dim - 1)
                    row++;
                if(row == z_dim)
                    row = 0;
            }
            else{
                dummy = tmp[j];
                // printf("index read by tmp is: %d, dma_length: %d \n", total_index, i+1);
            }
    	}
    }

    if(batch_total == 0)
    {
        //Copy X/P_pred into X/P in case we have more than one iteration
        row = 0;
        load_pred:for(unsigned i = 0; i < (x_dim + x_dim * x_dim); i++)
        {
#pragma HLS loop_tripcount max=42

            word_t tmp_pred;
            if(i < x_dim)
            {
                tmp_pred = X_pred[0][i];
                X[0][i] = tmp_pred;

// #ifndef __SYNTHESIS__
//                 printf("Copied value X_pred: %f \n", tmp);
// #endif

            }
            else
            {
                unsigned col = i - x_dim - row*x_dim;
                tmp_pred = P_pred[row][col];
                P[row][col] = tmp_pred;

                if(col == x_dim - 1)
                    row++;

// #ifndef __SYNTHESIS__
//                 printf("Copied value P_pred: %f \n", tmp);
// #endif
            }
        }
    }

}

void store(// word_t X_pred[1][X_MAX],
           // word_t P_pred[X_MAX][X_MAX],
           word_t _outbuff[SIZE_OUT_CHUNK_DATA],
           dma_word_t *out,
          /* <<--compute-params-->> */
           const unsigned iter,
           const unsigned x_dim,
           const unsigned z_dim,
	   dma_info_t &store_ctrl, int batch)
{
store_data:

    const unsigned length =
        round_up(x_dim + x_dim * x_dim, VALUES_PER_WORD);
    const unsigned store_offset =
        round_up(z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim , VALUES_PER_WORD) +
        round_up(z_dim, VALUES_PER_WORD) * (iter - 1);
    const unsigned out_offset = store_offset;
    const unsigned index = out_offset + length * batch;

    unsigned dma_length = length / VALUES_PER_WORD;
    unsigned dma_index = index / VALUES_PER_WORD;

#ifndef __SYNTHESIS__
    printf("dma_length %d\n", dma_length);
    printf("dma_index %d\n", dma_index);
#endif

    store_ctrl.index = dma_index;
    store_ctrl.length = dma_length;
    store_ctrl.size = SIZE_WORD_T;

    int row = 0;

    for (unsigned i = 0; i < dma_length; i++) {
#pragma HLS loop_tripcount max=21

    store_label1:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
	    out[dma_index + i].word[j] = _outbuff[i * VALUES_PER_WORD + j];

// #ifndef __SYNTHESIS__
//             int index = i * VALUES_PER_WORD + j;

//             if(index < x_dim)
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


void compute(word_t Z[1][Z_MAX],
             word_t X[1][X_MAX],
             word_t P[X_MAX][X_MAX],
             word_t F[X_MAX][X_MAX],
             word_t Q_kal[X_MAX][X_MAX],
             word_t R_kal[Z_MAX][Z_MAX],
             word_t H[Z_MAX][X_MAX],
             //word_t _inbuff[SIZE_IN_CHUNK_DATA],
             /* <<--compute-params-->> */
             const unsigned iter,
             const unsigned x_dim,
             const unsigned z_dim,
             word_t X_pred[1][X_MAX],
             word_t P_pred[X_MAX][X_MAX],
             word_t _outbuff[SIZE_OUT_CHUNK_DATA])
{
compute_data:

    // TODO implement compute functionality
    // const unsigned length =
    //     round_up(z_dim + x_dim + x_dim * x_dim * 3 + z_dim * z_dim + z_dim * x_dim, VALUES_PER_WORD) / 1;

    // for (int i = 0; i < length; i++)
    //     _outbuff[i] = _inbuff[i];

    // unsigned XX_size = x_dim*x_dim;
    // unsigned ZZ_size = z_dim*z_dim;
    // unsigned ZX_size = z_dim*x_dim;

    word_t inter1[X_MAX][X_MAX];
    word_t inter2[X_MAX][X_MAX];

    //Compute inter1 = F x P
    hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             X_MAX, X_MAX, X_MAX, X_MAX,
                             X_MAX, X_MAX,
                             TRAITS_A,
                             word_t, word_t> (F, P, inter1);

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
    // hls::print_matrix<X_MAX, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])inter2, "   ");
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

#ifndef __SYNTHESIS__
        // printf("inter2 = \n");
        // hls::print_matrix<X_MAX, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])inter2, "   ");
#endif

    //Compute S = H * inter2 * H.T + R
    word_t inter3[Z_MAX][X_MAX];
    word_t S[Z_MAX][Z_MAX];

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
            else if(i == j)//(i >= z_dim )
            {
                S[i][j] = 1.0;
            }
            else
                S[i][j] = 0.0;
        }

    //Compute S^-1 = S_inv - QR inverse
    word_t S_inv[Z_MAX][Z_MAX];
    int inv_ok = 0;

    //Initialize S_inv as Identity
    LOOP_S_inv_1:for(int i = 0; i < Z_MAX; i++)
    LOOP_S_inv_2:for(int j = 0; j < Z_MAX; j++)
        {
            if(i == j)//(i >= z_dim )
            {
                S_inv[i][j] = 1.0;
            }
            else
                S_inv[i][j] = 0.0;
        }

    //Use QR decomposition to compute inverse
    //hls::qr_inverse_top<Z_MAX,MY_CONFIG,word_t,word_t>(S,S_inv,inv_ok);
    //hls::qr_inverse<Z_MAX,word_t,word_t>(S,S_inv,inv_ok);
    inv_ok = inverse<Z_MAX>(S, S_inv);

#ifndef __SYNTHESIS__
    printf("inv_ok = %d\n", inv_ok);
    // printf("S_inv = \n");
    // hls::print_matrix<Z_MAX, Z_MAX, word_t, hls::NoTranspose>((word_t(*)[Z_MAX])S_inv, "   ");
#endif

    //Compute K = inter2 * H.T * S_inv
    word_t inter4[X_MAX][Z_MAX];
    word_t K[X_MAX][Z_MAX];

    //Compute inter4 = inter2 * H.T
    hls::matrix_multiply_top<hls::NoTranspose, hls::Transpose,
                             X_MAX, X_MAX, Z_MAX, X_MAX,
                             X_MAX, Z_MAX,
                             TRAITS_J,
                             word_t, word_t> (inter2, H, inter4);

#ifndef __SYNTHESIS__
    // printf("inter4 = \n");
    // hls::print_matrix<X_MAX, Z_MAX, word_t, hls::NoTranspose>((word_t(*)[Z_MAX])inter4, "   ");
#endif

    //Compute K = inter4 * S_inv
    hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             X_MAX, Z_MAX, Z_MAX, Z_MAX,
                             X_MAX, Z_MAX,
                             TRAITS_E,
                             word_t, word_t> (inter4, S_inv, K);

#ifndef __SYNTHESIS__
    // printf("K = \n");
    // hls::print_matrix<X_MAX, Z_MAX, word_t, hls::NoTranspose>((word_t(*)[Z_MAX])K, "   ");
#endif

    //Compute Y = Z - (H * F * X)
    //word_t inter5[Z_MAX][X_MAX];//Reuse inter3
    word_t Y[Z_MAX][1];

    //Compute inter5 = H * F
    hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             Z_MAX, X_MAX, X_MAX, X_MAX,
                             Z_MAX, X_MAX,
                             TRAITS_C,
                             word_t, word_t> (H, F, inter3);

    //Compute Y = inter5 * X
    hls::matrix_multiply_top<hls::NoTranspose, hls::Transpose,
                             Z_MAX, X_MAX, 1, X_MAX,
                             Z_MAX, 1,
                             TRAITS_F,
                             word_t, word_t> (inter3, X, Y);

    //Compute Y = Z - Y
    LOOP_Y:for(int i = 0; i < Z_MAX; i++)
    {
        if(i < z_dim)
        {
            word_t tmp = Y[i][0];
            Y[i][0] = Z[0][i] - tmp;
        }
        else
            Y[i][0] = 0.0;
    }

    //Compute final prediction state X_pred = F * X + K * Y
    //static word_t X_pred[1][X_MAX];//Set as PLMs from top
    word_t inter6[1][X_MAX];

    //Compute inter6 = K * Y
    hls::matrix_multiply_top<hls::Transpose, hls::Transpose,
                             Z_MAX, 1, X_MAX, Z_MAX,
                             1, X_MAX,
                             TRAITS_G,
                             word_t, word_t> (Y, K, inter6);

    //Compute X_pred = F * X (Compute X.T * F.T)
    hls::matrix_multiply_top<hls::NoTranspose, hls::Transpose,
                             1, X_MAX, X_MAX, X_MAX,
                             1, X_MAX,
                             TRAITS_H,
                             word_t, word_t> (X, F, X_pred);

    //Compute X_pred = X_pred + inter6
    LOOP_X_PRED:for(int i = 0; i < X_MAX; i++)
    {
        if(i < x_dim)
        {
            word_t tmp = X_pred[0][i];
            word_t tmp2 = inter6[0][i];
            X_pred[0][i] = tmp2 + tmp;
        }
        else
            X_pred[0][i] = 0.0;
    }

#ifndef __SYNTHESIS__
    // printf("X_pred = \n");
    // hls::print_matrix<1, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])X_pred, "   ");
#endif

    //Compute final prediction state covariance P_pred = (I - K * H) * inter2
    //static word_t P_pred[X_MAX][X_MAX];//Set as PLMs from top
    //word_t inter7[X_MAX][X_MAX];// Reuse inter1

    //Compute inter7 = K * H
    hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             X_MAX, Z_MAX, Z_MAX, X_MAX,
                             X_MAX, X_MAX,
                             TRAITS_I,
                             word_t, word_t> (K, H, inter1);

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
                word_t tmp = inter1[i][j];
                if(i == j)
                    inter1[i][j] = 1 - tmp;
                else
                    inter1[i][j] = 0 - tmp;
            }
            else
                inter1[i][j] = 0;
        }

    //Compute P_pred = inter7 * inter2
    hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose,
                             X_MAX, X_MAX, X_MAX, X_MAX,
                             X_MAX, X_MAX,
                             TRAITS_A,
                             word_t, word_t> (inter1, inter2, P_pred);

#ifndef __SYNTHESIS__
    // printf("P_pred = \n");
    // hls::print_matrix<X_MAX, X_MAX, word_t, hls::NoTranspose>((word_t(*)[X_MAX])P_pred, "   ");
#endif

    //Load output into outbuff
    int out_length = x_dim + x_dim * x_dim;
    int row = 0;

    LOOP_OUT:for(int i = 0; i < out_length; i++)
    {
#pragma HLS loop_tripcount max=42
        if(i < x_dim)
        {
            _outbuff[i] = X_pred[0][i];

#ifndef __SYNTHESIS__
            // printf("outbuff[%d] = %f\n", i, _outbuff[i]);
#endif
        }
        else
        {
            int col = i - x_dim - row*x_dim;
            _outbuff[i] = P_pred[row][col];

#ifndef __SYNTHESIS__
            // printf("outbuff[%d] = %.12f , P_pred[%d][%d] = %.12f \n", i, _outbuff[i], row, col, P_pred[row][col]);
#endif

            if(col == x_dim - 1)
                row++;

        }
    }

}


void top(dma_word_t *out, dma_word_t *in1,
         /* <<--params-->> */
	 const unsigned conf_info_iter,
	 const unsigned conf_info_x_dim,
	 const unsigned conf_info_z_dim,
	 dma_info_t &load_ctrl, dma_info_t &store_ctrl)
{

    /* <<--local-params-->> */
	 const unsigned iter = conf_info_iter;
	 const unsigned x_dim = conf_info_x_dim;
	 const unsigned z_dim = conf_info_z_dim;

         //static word_t _inbuff[SIZE_IN_CHUNK_DATA];

         static word_t Z[1][Z_MAX];
         static word_t X[1][X_MAX];
         static word_t P[X_MAX][X_MAX];
         static word_t F[X_MAX][X_MAX];
         static word_t Q_kal[X_MAX][X_MAX];
         static word_t R_kal[Z_MAX][Z_MAX];
         static word_t H[Z_MAX][X_MAX];

         static word_t _outbuff[SIZE_OUT_CHUNK_DATA];

         static word_t X_pred[1][X_MAX];
         static word_t P_pred[X_MAX][X_MAX];

    // Batching
batching:
    for (unsigned b = 0; b < iter; b++)
    {
#pragma HLS loop_tripcount max=1
        // Chunking
    go:
        // for (int c = 0; c < 1; c++)
        // {

        load(Z,
             X,
             P,
             F,
             Q_kal,
             R_kal,
             H,
             X_pred,
             P_pred,
             //_inbuff,
             in1,
             /* <<--args-->> */
             iter,
             x_dim,
             z_dim,
             load_ctrl, b);

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
                iter,
                x_dim,
                z_dim,
                X_pred,
                P_pred,
                _outbuff);
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
              iter,
              x_dim,
              z_dim,
              store_ctrl, b);

        // }
    }
}


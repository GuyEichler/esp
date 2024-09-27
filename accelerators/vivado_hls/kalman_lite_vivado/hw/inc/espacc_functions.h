
#include "espacc.h"

//helper functions
template<unsigned DIM, typename T>
int inverse(T new_mat[DIM][DIM], T out[DIM][DIM], unsigned real_dim)
/* void inverse(T new_mat[DIM][DIM], T out[DIM][DIM], unsigned real_dim, int& inv_ok) */
{
//#pragma HLS inline
inverse_func:

	//int i,j,k;

    if(DIM == 2){
        T a = new_mat[0][0];
        T b = new_mat[0][1];
        T c = new_mat[1][0];
        T d = new_mat[1][1];

        T det = (a * d) - (b * c);

        if (det == 0) {
#ifndef __SYNTHESIS__
            printf("The matrix is not invertible.\n");
#endif
            return 1;
	    /* inv_ok = 1; */
            /* return; */
        }

        out[0][0] = d / det;
        out[0][1] = (-1) * b / det;
        out[1][0] = (-1) * c / det;
        out[1][1] = a / det;
    }
    else
    {
        /* Applying Gauss Jordan Elimination */
        /* main_loop:for(i = 0; i < DIM; i++) */
        main_loop:for(int i = 0; i < real_dim; i++)
        {
         #pragma HLS loop_tripcount max=164
	 /* #pragma HLS dataflow */
	    /* loop_1:for(j = 0; j < DIM; j++) */
	    loop_1:for(int j = 0; j < real_dim; j++)
            {
            #pragma HLS loop_tripcount max=164
            /* #pragma HLS dataflow */
                if(i != j)
                {
                    T ratio;
		    T new_ji = new_mat[j][i];
		    T new_ii = new_mat[i][i];
		    /* T new_ii_inv = 1 / new_ii; */
                    ratio = new_ji / new_ii;
                    /* ratio = new_ji * new_ii_inv; */

		    T new_jj = 0.0; //new_mat[j][j];
		    T new_jj_inv = 0.0;

		    /* if(i == DIM - 1) */
		    if(i == real_dim - 1)
		    {
			new_jj = new_mat[j][j];
		        T new_ij = new_mat[i][j];
		        T ratio_new_ij = ratio * new_ij;
			T new_val = new_jj - ratio_new_ij;
		        new_mat[j][j] = new_val;
			new_jj = new_val;
			new_jj_inv = 1 / new_val;

		    }


		    /* loop_2:for(k = 0; k < DIM; k++) */
		    loop_2:for(int k = 0; k < real_dim; k++)
                    {
                    #pragma HLS loop_tripcount max=164
                    /* #pragma HLS dataflow */
                    #pragma HLS PIPELINE II=2
                    #pragma HLS dependence variable=out inter false
                    #pragma HLS dependence variable=new_mat inter false
                    /* #pragma HLS array_partition \ */
		    /* 	    	variable=new_mat cyclic factor=2 dim=1 */
			T new_jk = new_mat[j][k];
			T new_ik = new_mat[i][k];
			T new_jk_val = new_jk - ratio*new_ik;

			T out_jk = out[j][k];
			T out_ik = out[i][k];

			T out_jk_val = out_jk - ratio*out_ik;

                        /* if(i == DIM-1){ */
                        if(i == real_dim - 1){
                            if(k != j || k == 0){
				new_mat[j][k] = new_jk_val;
                            }

			    out[j][k] = out_jk_val * new_jj_inv;

                        }
                        else{

                            new_mat[j][k] = new_jk_val;

			    out[j][k] = out_jk_val;

		        }
		    }
		}
            }
        }



	/* Row Operation to Make Principal Diagonal to 1 */
        T new_last = new_mat[real_dim-1][real_dim-1];
        /* loop_3:for(i = 0; i < DIM; i++) */
        loop_3:for(int i = 0; i < real_dim; i++)
        {
        #pragma HLS loop_tripcount max=164
        #pragma HLS PIPELINE
        #pragma HLS dependence variable=out inter false
	    /* T out_i = out[DIM-1][i]; */
	    /* T new_last = new_mat[DIM-1][DIM-1]; */
            /* out[DIM-1][i] = out_i / new_last; */
	    T out_i = out[real_dim-1][i];
	    /* T new_last = new_mat[real_dim-1][real_dim-1]; */
            out[real_dim-1][i] = out_i / new_last;
        }

    }

    return 0;
    /* inv_ok = 0; */
    /* return; */
}

template<unsigned DIM, typename T>
int inverse2(T new_mat[DIM][DIM], T out[DIM][DIM], unsigned real_dim)
/* void inverse(T new_mat[DIM][DIM], T out[DIM][DIM], unsigned real_dim, int& inv_ok) */
{
//#pragma HLS inline
inverse_func:

	//int i,j,k;

    if(DIM == 2){
        T a = new_mat[0][0];
        T b = new_mat[0][1];
        T c = new_mat[1][0];
        T d = new_mat[1][1];

        T det = (a * d) - (b * c);

        if (det == 0) {
#ifndef __SYNTHESIS__
            printf("The matrix is not invertible.\n");
#endif
            return 1;
	    /* inv_ok = 1; */
            /* return; */
        }

        out[0][0] = d / det;
        out[0][1] = (-1) * b / det;
        out[1][0] = (-1) * c / det;
        out[1][1] = a / det;
    }
    else
    {
        /* Applying Gauss Jordan Elimination */
        /* main_loop:for(i = 0; i < DIM; i++) */
        main_loop:for(int i = 0; i < real_dim-1; i++)
        {
         #pragma HLS loop_tripcount max=164
	 /* #pragma HLS dataflow */
	    /* loop_1:for(j = 0; j < DIM; j++) */
	    loop_1:for(int j = 0; j < real_dim; j++)
            {
            #pragma HLS loop_tripcount max=164
            /* #pragma HLS dataflow */
                if(i != j)
                {
                    T ratio;
		    T new_ji = new_mat[j][i];
		    T new_ii = new_mat[i][i];
		    /* T new_ii_inv = 1 / new_ii; */
                    ratio = new_ji / new_ii;
                    /* ratio = new_ji * new_ii_inv; */

		    T new_jj = 0.0; //new_mat[j][j];
		    T new_jj_inv = 0.0;

		    /* if(i == DIM - 1) */
		    /* if(i == real_dim - 1) */
		    /* { */
		    /* 	new_jj = new_mat[j][j]; */
		    /*     T new_ij = new_mat[i][j]; */
		    /*     T ratio_new_ij = ratio * new_ij; */
		    /* 	T new_val = new_jj - ratio_new_ij; */
		    /*     new_mat[j][j] = new_val; */
		    /* 	new_jj = new_val; */
		    /* 	new_jj_inv = 1 / new_val; */

		    /* } */


		    /* loop_2:for(k = 0; k < DIM; k++) */
		    loop_2:for(int k = 0; k < real_dim; k++)
                    {
                    #pragma HLS loop_tripcount max=164
                    /* #pragma HLS dataflow */
                    #pragma HLS PIPELINE II=1

			T new_jk = new_mat[j][k];
			T new_ik = new_mat[i][k];
			T new_jk_val = new_jk - ratio*new_ik;

			T out_jk = out[j][k];
			T out_ik = out[i][k];

			T out_jk_val = out_jk - ratio*out_ik;

                        /* if(i == DIM-1){ */
                        /* if(i == real_dim - 1){ */
                        /*     if(k != j || k == 0){ */
			/* 	new_mat[j][k] = new_jk_val; */
                        /*     } */

			/*     out[j][k] = out_jk_val * new_jj_inv; */

                        /* } */
                        /* else{ */

                            new_mat[j][k] = new_jk_val;

			    out[j][k] = out_jk_val;

		        }
		    }
		}
            }
        /* } */

        /* main_loop:for(int i = 0; i < real_dim; i++) */
        /* { */
        /*  #pragma HLS loop_tripcount max=164 */
	 /* #pragma HLS dataflow */
	    /* loop_1:for(j = 0; j < DIM; j++) */
	    loop_4:for(int j = 0; j < real_dim; j++)
            {
            #pragma HLS loop_tripcount max=164
            /* #pragma HLS dataflow */
		int i = real_dim - 1;
                if(i != j)
                {
                    T ratio;
		    T new_ji = new_mat[j][i];
		    T new_ii = new_mat[i][i];
		    /* T new_ii_inv = 1 / new_ii; */
                    ratio = new_ji / new_ii;
                    /* ratio = new_ji * new_ii_inv; */

		    T new_jj = 0.0; //new_mat[j][j];
		    T new_jj_inv = 0.0;

		    /* if(i == DIM - 1) */
		    /* if(i == real_dim - 1) */
		    /* { */
			new_jj = new_mat[j][j];
		        T new_ij = new_mat[i][j];
		        T ratio_new_ij = ratio * new_ij;
			T new_val = new_jj - ratio_new_ij;
		        new_mat[j][j] = new_val;
			new_jj = new_val;
			new_jj_inv = 1 / new_val;

		    /* } */


		    /* loop_2:for(k = 0; k < DIM; k++) */
		    loop_5:for(int k = 0; k < j; k++)
                    {
                    #pragma HLS loop_tripcount max=164
                    /* #pragma HLS dataflow */
                    #pragma HLS PIPELINE II=1

			T new_jk = new_mat[j][k];
			T new_ik = new_mat[i][k];
			T new_jk_val = new_jk - ratio*new_ik;

			T out_jk = out[j][k];
			T out_ik = out[i][k];

			T out_jk_val = out_jk - ratio*out_ik;

                        /* if(i == DIM-1){ */
                        /* if(i == real_dim - 1){ */
                            /* if(k != j || k == 0){ */
				new_mat[j][k] = new_jk_val;
                            /* } */

			    out[j][k] = out_jk_val * new_jj_inv;

                        /* } */
                        /* else{ */

                        /*     new_mat[j][k] = new_jk_val; */

			/*     out[j][k] = out_jk_val; */

		        /* } */
		    }


                    //k == j
			T new_jk_eq = new_mat[j][j];
			T new_ik_eq = new_mat[i][j];
			T new_jk_val_eq = new_jk_eq - ratio*new_ik_eq;

			T out_jk_eq = out[j][j];
			T out_ik_eq = out[i][j];

			T out_jk_val_eq = out_jk_eq - ratio*out_ik_eq;

                        /* if(i == DIM-1){ */
                        /* if(i == real_dim - 1){ */
                            /* if(k != j || k == 0){ */
				new_mat[j][j] = new_jk_val_eq;
                            /* } */

			    out[j][j] = out_jk_val_eq * new_jj_inv;


		    loop_6:for(int k = j+1; k < real_dim; k++)
                    {
                    #pragma HLS loop_tripcount max=164
                    /* #pragma HLS dataflow */
                    #pragma HLS PIPELINE II=1

			T new_jk = new_mat[j][k];
			T new_ik = new_mat[i][k];
			T new_jk_val = new_jk - ratio*new_ik;

			T out_jk = out[j][k];
			T out_ik = out[i][k];

			T out_jk_val = out_jk - ratio*out_ik;

                        /* if(i == DIM-1){ */
                        /* if(i == real_dim - 1){ */
                            /* if(k != j || k == 0){ */
				new_mat[j][k] = new_jk_val;
                            /* } */

			    out[j][k] = out_jk_val * new_jj_inv;

                        /* } */
                        /* else{ */

                        /*     new_mat[j][k] = new_jk_val; */

			/*     out[j][k] = out_jk_val; */

		        /* } */
		    }





	       }
	}
        /* } */



	/* Row Operation to Make Principal Diagonal to 1 */
        /* loop_3:for(i = 0; i < DIM; i++) */
        loop_3:for(int i = 0; i < real_dim; i++)
        {
        #pragma HLS loop_tripcount max=164
        #pragma HLS PIPELINE
	    /* T out_i = out[DIM-1][i]; */
	    /* T new_last = new_mat[DIM-1][DIM-1]; */
            /* out[DIM-1][i] = out_i / new_last; */
	    T out_i = out[real_dim-1][i];
	    T new_last = new_mat[real_dim-1][real_dim-1];
            out[real_dim-1][i] = out_i / new_last;
        }

	}

    return 0;
    /* inv_ok = 0; */
    /* return; */
}

template<unsigned DIM, typename T>
int inverse_partial(T new_mat[DIM][DIM], T out[DIM][DIM], unsigned real_dim)
/* void inverse(T new_mat[DIM][DIM], T out[DIM][DIM], unsigned real_dim, int& inv_ok) */
{
//#pragma HLS inline
inverse_func:

	//int i,j,k;

    if(DIM == 2){
        T a = new_mat[0][0];
        T b = new_mat[0][1];
        T c = new_mat[1][0];
        T d = new_mat[1][1];

        T det = (a * d) - (b * c);

        if (det == 0) {
#ifndef __SYNTHESIS__
            printf("The matrix is not invertible.\n");
#endif
            return 1;
	    /* inv_ok = 1; */
            /* return; */
        }

        out[0][0] = d / det;
        out[0][1] = (-1) * b / det;
        out[1][0] = (-1) * c / det;
        out[1][1] = a / det;
    }
    else
    {
        /* Applying Gauss Jordan Elimination */
        /* main_loop:for(i = 0; i < DIM; i++) */
        main_loop:for(int i = 0; i < real_dim/2; i++)
        {
         #pragma HLS loop_tripcount max=164/2
	 /* #pragma HLS dataflow */
	    /* loop_1:for(j = 0; j < DIM; j++) */
	    loop_1:for(int j = 0; j < real_dim; j++)
            {
            #pragma HLS loop_tripcount max=164
            /* #pragma HLS dataflow */
                if(i != j)
                {
                    T ratio;
		    T new_ji = new_mat[j][i];
		    T new_ii = new_mat[i][i];
		    /* T new_ii_inv = 1 / new_ii; */
                    ratio = new_ji / new_ii;
                    /* ratio = new_ji * new_ii_inv; */

		    T new_jj = 0.0; //new_mat[j][j];
		    T new_jj_inv = 0.0;

		    /* if(i == DIM - 1) */
		    if(i == real_dim - 1)
		    {
			new_jj = new_mat[j][j];
		        T new_ij = new_mat[i][j];
		        T ratio_new_ij = ratio * new_ij;
			T new_val = new_jj - ratio_new_ij;
		        new_mat[j][j] = new_val;
			new_jj = new_val;
			new_jj_inv = 1 / new_val;

		    }


		    /* loop_2:for(k = 0; k < DIM; k++) */
		    loop_2:for(int k = 0; k < real_dim; k++)
                    {
                    #pragma HLS loop_tripcount max=164
                    /* #pragma HLS dataflow */
                    #pragma HLS PIPELINE II=2

			T new_jk = new_mat[j][k];
			T new_ik = new_mat[i][k];
			T new_jk_val = new_jk - ratio*new_ik;

			T out_jk = out[j][k];
			T out_ik = out[i][k];

			T out_jk_val = out_jk - ratio*out_ik;

                        /* if(i == DIM-1){ */
                        if(i == real_dim - 1){
                            if(k != j || k == 0){
				new_mat[j][k] = new_jk_val;
                            }

			    out[j][k] = out_jk_val * new_jj_inv;

                        }
                        else{

                            new_mat[j][k] = new_jk_val;

			    out[j][k] = out_jk_val;

		        }
		    }
		}
            }
        }



	/* Row Operation to Make Principal Diagonal to 1 */
        /* loop_3:for(i = 0; i < DIM; i++) */
        loop_3:for(int i = 0; i < real_dim; i++)
        {
        #pragma HLS loop_tripcount max=164
        #pragma HLS PIPELINE
	    /* T out_i = out[DIM-1][i]; */
	    /* T new_last = new_mat[DIM-1][DIM-1]; */
            /* out[DIM-1][i] = out_i / new_last; */
	    T out_i = out[real_dim-1][i];
	    T new_last = new_mat[real_dim-1][real_dim-1];
            out[real_dim-1][i] = out_i / new_last;
        }

    }

    return 0;
    /* inv_ok = 0; */
    /* return; */
}


void compute_int2(word_t inter2[X_MAX][X_MAX], word_t Q_kal[X_MAX][X_MAX], unsigned x_dim)
{
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
}

void compute_sr(word_t S[Z_MAX][Z_MAX], word_t R_kal[Z_MAX][Z_MAX], unsigned z_dim)
{
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
}

void init_s_inv(word_t S_inv[Z_MAX][Z_MAX], unsigned z_dim)
{
    LOOP_S_inv_1:for(int i = 0; i < Z_MAX; i++)
    LOOP_S_inv_2:for(int j = 0; j < Z_MAX; j++)
        {
            if(i == j && i < z_dim)
            {
                S_inv[i][j] = 1.0;
            }
            else
                S_inv[i][j] = 0.0;
        }
}

void compute_y(word_t Y[Z_MAX][1], word_t Z[CHUNK_MAX][Z_MAX], unsigned z_dim, unsigned curr_chunk)
{
    LOOP_Y:for(int i = 0; i < Z_MAX; i++)
    {
        if(i < z_dim)
        {
            word_t tmp = Y[i][0];
            Y[i][0] = Z[curr_chunk][i] - tmp;
// #ifndef __SYNTHESIS__
//                 printf("Computed with Z[%d][%d]: %f\n", curr_chunk, i, Z[curr_chunk][i]);
//                 printf("Y[%d][%d]: %f\n", i, 0, Y[i][0]);
// #endif
        }
        else
            Y[i][0] = 0.0;
    }
}

void compute_x_pred(word_t X_pred[1][X_MAX], word_t inter6[1][X_MAX], unsigned x_dim)
{
    LOOP_X_PRED:for(int i = 0; i < X_MAX; i++)
    {
        if(i < x_dim)
        {
            word_t tmp = X_pred[0][i];
            word_t tmp2 = inter6[0][i];
            X_pred[0][i] = tmp2 + tmp;
// #ifndef __SYNTHESIS__
//                 printf("Compute value X_pred: %f \n", X_pred[0][i]);
// #endif
        }
        else
            X_pred[0][i] = 0.0;
    }
}

void compute_int7(word_t inter7[X_MAX][X_MAX], unsigned x_dim)
{
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
}

void compute_out(word_t _outbuff[SIZE_OUT_CHUNK_DATA], word_t X_pred[1][X_MAX], word_t P_pred[X_MAX][X_MAX], unsigned x_dim, unsigned curr_chunk)
{

    int out_length = x_dim + x_dim * x_dim;
    int row = 0;
    LOOP_OUT:for(int i = 0; i < X_MAX + X_MAX * X_MAX; i++)
    {
//#pragma HLS loop_tripcount max=42
        if(i < x_dim)
        {
            _outbuff[i + curr_chunk * out_length] = X_pred[0][i];

// #ifndef __SYNTHESIS__
//             printf("compute outbuff[%d] = %f\n", i, _outbuff[i + curr_chunk * out_length]);
// #endif
        }
        else if(i < x_dim + x_dim * x_dim)
        {
            int col = i - x_dim - row*x_dim;
            _outbuff[i + curr_chunk * out_length] = P_pred[row][col];

#ifndef __SYNTHESIS__
            // printf("outbuff[%d] = %.12f , P_pred[%d][%d] = %.12f \n", i, _outbuff[i], row, col, P_pred[row][col]);
#endif

            if(col == x_dim - 1)
                row++;

        }
    }
}

template<unsigned X_DIM, unsigned Z_DIM>
void mat_mul(word_t A[X_DIM][Z_DIM], word_t B[Z_DIM][Z_DIM], word_t C[X_DIM][Z_DIM], unsigned x_dim, unsigned z_dim)
{

LOOP_K1:for (int i = 0; i < x_dim; i++) {
         #pragma HLS loop_tripcount max=164
    /* #pragma HLS unroll factor=1 */
			    /* K[i][j] = 0; */
    LOOP_K3:for (int k = 0; k < z_dim; k++) {
    /* #pragma HLS unroll factor=1 */
         #pragma HLS loop_tripcount max=164
	    word_t A_ik = A[i][k];
	    LOOP_K2:for (int j = 0; j < z_dim; j++) {
                            #pragma HLS loop_tripcount max=164
                            #pragma HLS array_partition \
			    	variable=C cyclic factor=8 dim=2
                            /* #pragma HLS array_partition \ */
			    /* 	variable=A cyclic factor=4 dim=2 */
                            #pragma HLS array_partition \
			    	variable=B cyclic factor=8 dim=2
                            #pragma HLS dependence variable=C inter false
                            #pragma HLS dependence variable=B inter false
                            #pragma HLS PIPELINE II=1
                            #pragma HLS unroll factor=8
				    /* if(i < x_dim && j < z_dim && k < z_dim) */
					    if(k == 0)
						    /* C[i][j] = A[i][k] * B[k][j]; */
						    C[i][j] = A_ik * B[k][j];
					    else{
						    /* C[i][j] += A[i][k] * B[k][j]; */
						    C[i][j] += A_ik * B[k][j];
					    }
	    }
		}
	}

}

template<unsigned X_DIM, unsigned Z_DIM>
void mat_mul_sub(word_t A[X_DIM][Z_DIM], word_t B[Z_DIM][Z_DIM], word_t C[X_DIM][Z_DIM], unsigned x_dim, unsigned z_dim, word_t sub)
{

LOOP_K1:for (int i = 0; i < X_DIM; i++) {
    /* #pragma HLS unroll factor=1 */
			    /* K[i][j] = 0; */
    LOOP_K3:for (int k = 0; k < Z_DIM; k++) {
    /* #pragma HLS unroll factor=1 */
	    word_t A_ik = A[i][k];
	    LOOP_K2:for (int j = 0; j < Z_DIM; j++) {
                            #pragma HLS array_partition \
			    	variable=C cyclic factor=8 dim=2
                            /* #pragma HLS array_partition \ */
			    /* 	variable=A cyclic factor=4 dim=2 */
                            #pragma HLS array_partition \
			    	variable=B cyclic factor=8 dim=2
                            #pragma HLS dependence variable=C inter false
                            #pragma HLS dependence variable=B inter false
                            #pragma HLS PIPELINE II=1
                            #pragma HLS unroll factor=8
				    if(i < x_dim && j < z_dim && k < z_dim)
					    if(k == 0)
						    /* C[i][j] = A[i][k] * B[k][j]; */
						    C[i][j] = A_ik * B[k][j];
					    else{
						    /* C[i][j] += A[i][k] * B[k][j]; */
						    C[i][j] += A_ik * B[k][j];
					    }

		if(k == z_dim - 1)
		    if(i == j)
			    C[i][j] = sub - C[i][j];
		    else
			    C[i][j] = 0 - C[i][j];
	    }

	    }
    }

}

template<unsigned DIM>
void iterative_inverse(word_t new_mat[DIM][DIM], word_t out[DIM][DIM], word_t new_out_final[DIM][DIM], unsigned z_dim)
{
	word_t new_out[DIM][DIM];// = {{0}};
	//float new_out_final[DIM][DIM] = {{0}};
	//float I2[DIM][DIM];

	//A*X_n
	mat_mul<DIM, DIM>(new_mat, out, new_out, z_dim, z_dim);
	/* hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose, */
        /*                      Z_MAX, Z_MAX, Z_MAX, Z_MAX, */
        /*                      Z_MAX, Z_MAX, */
        /*                      TRAITS_K, */
        /*                      word_t, word_t> (new_mat, out, new_out); */

	//2I - A*X_n
	for(int i = 0; i < z_dim; i++)
         #pragma HLS loop_tripcount max=164
		for(int j = 0; j < z_dim; j++){
         #pragma HLS loop_tripcount max=164
                    #pragma HLS PIPELINE
                    #pragma HLS dependence variable=new_out inter false
			/* if(i < z_dim && j < z_dim) */
				if(i == j){
					new_out[i][j] = 2 - new_out[i][j];
				}
				else{
					new_out[i][j] = 0 - new_out[i][j];
				}
	}

	//X_n*(2I-A*X_n)
	mat_mul<DIM>(out, new_out, new_out_final, z_dim, z_dim);
	/* hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose, */
        /*                      Z_MAX, Z_MAX, Z_MAX, Z_MAX, */
        /*                      Z_MAX, Z_MAX, */
        /*                      TRAITS_K, */
        /*                      word_t, word_t> (out, new_out, new_out_final); */
    /* 	printf("final matrix = \n"); */
    /* hls::print_matrix<Z_MAX, Z_MAX, word_t, hls::NoTranspose>((word_t(*)[Z_MAX])new_out_final, "   "); */
}

template<unsigned DIM>
void iterative_inverse_lite(word_t new_mat[DIM][DIM], word_t out[DIM][DIM], unsigned z_dim)
{
	word_t new_out[DIM][DIM];// = {{0}};
	//float new_out_final[DIM][DIM] = {{0}};
	//float I2[DIM][DIM];

	//A*X_n
	mat_mul<DIM, DIM>(new_mat, out, new_out, z_dim, z_dim);
	/* hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose, */
        /*                      Z_MAX, Z_MAX, Z_MAX, Z_MAX, */
        /*                      Z_MAX, Z_MAX, */
        /*                      TRAITS_K, */
        /*                      word_t, word_t> (new_mat, out, new_out); */

	//2I - A*X_n
	for(int i = 0; i < z_dim; i++)
         #pragma HLS loop_tripcount max=164
		for(int j = 0; j < z_dim; j++){
         #pragma HLS loop_tripcount max=164
                    #pragma HLS PIPELINE
                    #pragma HLS dependence variable=new_out inter false
                    #pragma HLS dependence variable=new_mat inter false
                    #pragma HLS dependence variable=out inter false
                    #pragma HLS unroll factor=8
                    #pragma HLS array_partition variable=new_mat cyclic factor=8 dim=2
			/* if(i < z_dim && j < z_dim) */
   			        word_t tmp2 = out[i][j];
   			        word_t tmp = new_out[i][j];
				word_t result;
				if(i == j){
					result = 2 - tmp;
				}
				else{
					result = 0 - tmp;
				}
				new_out[i][j] = result;
				new_mat[i][j] = tmp2;//Reuse new_mat
	}

	//X_n*(2I-A*X_n)
	mat_mul<DIM>(new_mat, new_out, out, z_dim, z_dim);
	/* hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose, */
        /*                      Z_MAX, Z_MAX, Z_MAX, Z_MAX, */
        /*                      Z_MAX, Z_MAX, */
        /*                      TRAITS_K, */
        /*                      word_t, word_t> (out, new_out, new_out_final); */
    /* 	printf("final matrix = \n"); */
    /* hls::print_matrix<Z_MAX, Z_MAX, word_t, hls::NoTranspose>((word_t(*)[Z_MAX])new_out_final, "   "); */
}

template<unsigned DIM, unsigned ITER>
	void many_iterative_inverse(word_t new_mat[DIM][DIM], word_t out[DIM][DIM], unsigned z_dim)
{

	word_t new_out_final[DIM][DIM];
	for(int k = 0; k < ITER; k++){

		//std::cout << "ITERATION " << k << std::endl;

		//float next_out[DIM][DIM] = {{0}};
		iterative_inverse<DIM>(new_mat, out, new_out_final, z_dim);
		//iterative_inverse<N>(mat, rand_inv, next_out);

		//Initialize inverse as the new iteration
		for(int i = 0; i < DIM; i++)
			for(int j = 0; j < DIM; j++)
                        #pragma HLS PIPELINE
				out[i][j] = new_out_final[i][j];
		//rand_inv[i][j] = next_out[i][j];

	}
}

template<unsigned DIM>
void iterative_inverse2(word_t new_mat[DIM][DIM], word_t out[DIM][DIM], word_t final_out[DIM][DIM],word_t final_out2[DIM][DIM], unsigned z_dim)
{
/* #pragma HLS dataflow */
    iterative_inverse<DIM>(new_mat, out, final_out, z_dim);
    iterative_inverse<DIM>(new_mat, final_out, final_out2, z_dim);

}

template<unsigned DIM>
void iterative_inverse3(word_t new_mat[DIM][DIM], word_t out[DIM][DIM], word_t final_out[DIM][DIM],word_t final_out2[DIM][DIM], unsigned z_dim)
{
/* #pragma HLS dataflow */
    iterative_inverse<DIM>(new_mat, out, final_out, z_dim);
    iterative_inverse<DIM>(new_mat, final_out, final_out2, z_dim);
    iterative_inverse<DIM>(new_mat, final_out2, final_out, z_dim);
}

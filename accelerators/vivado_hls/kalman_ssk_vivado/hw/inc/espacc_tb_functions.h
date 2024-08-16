
//#include "espacc.h"

//helper functions
template<unsigned DIM, typename T>
int inverse_tb(T new_mat[DIM][DIM], T out[DIM][DIM], unsigned real_dim)
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

template<unsigned X_DIM, unsigned Z_DIM>
void mat_mul_tb(word_t A[X_DIM][Z_DIM], word_t B[Z_DIM][Z_DIM], word_t C[X_DIM][Z_DIM], unsigned x_dim, unsigned z_dim)
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
void mat_mul_sub_tb(word_t A[X_DIM][Z_DIM], word_t B[Z_DIM][Z_DIM], word_t C[X_DIM][Z_DIM], unsigned x_dim, unsigned z_dim, word_t sub)
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
void iterative_inverse_tb(word_t new_mat[DIM][DIM], word_t out[DIM][DIM], word_t new_out_final[DIM][DIM], unsigned z_dim)
{
	word_t new_out[DIM][DIM];// = {{0}};
	//float new_out_final[DIM][DIM] = {{0}};
	//float I2[DIM][DIM];

	//A*X_n
	mat_mul_tb<DIM, DIM>(new_mat, out, new_out, z_dim, z_dim);
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
	mat_mul_tb<DIM>(out, new_out, new_out_final, z_dim, z_dim);
	/* hls::matrix_multiply_top<hls::NoTranspose, hls::NoTranspose, */
        /*                      Z_MAX, Z_MAX, Z_MAX, Z_MAX, */
        /*                      Z_MAX, Z_MAX, */
        /*                      TRAITS_K, */
        /*                      word_t, word_t> (out, new_out, new_out_final); */
    /* 	printf("final matrix = \n"); */
    /* hls::print_matrix<Z_MAX, Z_MAX, word_t, hls::NoTranspose>((word_t(*)[Z_MAX])new_out_final, "   "); */
}


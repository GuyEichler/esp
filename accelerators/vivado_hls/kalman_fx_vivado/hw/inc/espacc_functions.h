
#include "espacc.h"

//helper functions
template<unsigned DIM, typename T>
int inverse(T new_mat[DIM][DIM], T out[DIM][DIM], unsigned real_dim)
{
//#pragma HLS inline
inverse_func:

    int i,j,k;

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
        main_loop:for(i = 0; i < real_dim; i++)
        {
         #pragma HLS loop_tripcount max=164
	    /* loop_1:for(j = 0; j < DIM; j++) */
	    loop_1:for(j = 0; j < real_dim; j++)
            {
            #pragma HLS loop_tripcount max=164
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
		    loop_2:for(k = 0; k < real_dim; k++)
                    {
                    #pragma HLS loop_tripcount max=164
#pragma HLS dependence variable=new_mat inter false
#pragma HLS dependence variable=out inter false
/* #pragma HLS array_partition variable=new_mat cyclic factor=32 dim=2 */
/* #pragma HLS array_partition variable=out cyclic factor=32 dim=2 */
/* #pragma HLS array_partition variable=new_mat cyclic factor=16 dim=1 */
/* #pragma HLS array_partition variable=out cyclic factor=16 dim=1 */
                    #pragma HLS PIPELINE II=2

			T new_jk = new_mat[j][k];
			T new_ik = new_mat[i][k];
			T ratio_new_ik = ratio*new_ik;
			T new_jk_val = new_jk - ratio_new_ik;

			T out_jk = out[j][k];
			T out_ik = out[i][k];
			T ratio_out_ik = ratio*out_ik;

			T out_jk_val = out_jk - ratio_out_ik;

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
        loop_3:for(i = 0; i < real_dim; i++)
        {
        #pragma HLS loop_tripcount max=164
#pragma HLS dependence variable=out inter false
        #pragma HLS PIPELINE
	    /* T out_i = out[DIM-1][i]; */
	    /* T new_last = new_mat[DIM-1][DIM-1]; */
            /* out[DIM-1][i] = out_i / new_last; */
	    T out_i = out[real_dim-1][i];
	    /* T new_last = new_mat[real_dim-1][real_dim-1]; */
            out[real_dim-1][i] = out_i / new_last;
        }

    }

    return 0;
}

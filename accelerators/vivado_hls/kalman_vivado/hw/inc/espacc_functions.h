
#include "espacc.h"

//helper functions
template<unsigned DIM>
int inverse(word_t new_mat[DIM][DIM], word_t out[DIM][DIM])
{
inverse_func:

    int i,j,k;

    if(DIM == 2){
        word_t a = new_mat[0][0];
        word_t b = new_mat[0][1];
        word_t c = new_mat[1][0];
        word_t d = new_mat[1][1];

        word_t det = (a * d) - (b * c);

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
        main_loop:for(i = 0; i < DIM; i++)
        {
	    loop_1:for(j = 0; j < DIM; j++)
            {
                if(i != j)
                {
                    word_t ratio;
		    word_t new_ji = new_mat[j][i];
		    word_t new_ii = new_mat[i][i];
		    /* word_t new_ii_inv = 1 / new_ii; */
                    ratio = new_ji / new_ii;
                    /* ratio = new_ji * new_ii_inv; */

		    word_t new_jj = new_mat[j][j];
		    word_t new_jj_inv;
		    if(i == DIM - 1)
		    {
		        word_t new_ij = new_mat[i][j];
		        word_t ratio_new_ij = ratio * new_ij;
			word_t new_val = new_jj - ratio_new_ij;
		        new_mat[j][j] = new_val;
			new_jj = new_val;
			new_jj_inv = 1 / new_val;

			/* word_t new_j0 = new_mat[j][0]; */
			/* word_t new_i0 = new_mat[i][0]; */
			/* word_t ratio_new_i0 = ratio * new_i0; */
			/* //new_j0 = new_j0 - ratio_new_i0; */
			/* new_mat[j][0] = new_j0 - ratio*new_i0; */

		    }
		    //new_jj = new_jj - ratio_new_jj;



		    loop_2:for(k = 0; k < DIM; k++)
                    {
                    #pragma HLS PIPELINE II=2

			word_t new_jk = new_mat[j][k];
			word_t new_ik = new_mat[i][k];
			word_t new_jk_val = new_jk - ratio*new_ik;

			word_t out_jk = out[j][k];
			word_t out_ik = out[i][k];
			/* out[j][k] = out_jk - ratio*out_ik; */
			word_t out_jk_val = out_jk - ratio*out_ik;

                        if(i == DIM-1){
                            /* if(k == 0){//Calc the diagonal element first */
                            /*     /\* word_t new_jj = new_mat[j][j]; *\/ */
                            /*     //new_mat[j][j] = new_jj - ratio_new_ij; */
                            /*     //out[j][j] = (out[j][j] - ratio*out[i][j]) / */
                            /*     //new_mat[j][j]; */
                            /* } */
                            /* else if(k == j){ */
			    /*     /\* word_t new_j0 = new_mat[j][0]; *\/ */
			    /*     /\* word_t new_i0 = new_mat[i][0]; *\/ */
                            /*     //new_mat[j][0] = new_j0 - ratio*new_i0; */
                            /*     //out[j][0] = (out[j][0] - ratio*out[i][0]) / */
                            /*     //new_mat[j][j]; */
                            /* } */
                            if(k != j || k == 0){
			        /* word_t new_jk = new_mat[j][k]; */
			        /* word_t new_ik = new_mat[i][k]; */
                                /* new_mat[j][k] = new_jk - ratio*new_ik; */
				new_mat[j][k] = new_jk_val;
                                //out[j][k] = (out[j][k] - ratio*out[i][k]) /
                                //new_mat[j][j];
                            }

			    //word_t out_jk = out[j][k];
			    /* word_t out_ik = out[i][k]; */
			    /* //word_t local_jj = new_mat[j][j]; */
			    /* //new_jj = new_mat[j][j]; */
                            /* /\* out[j][k] = (out_jk - ratio*out_ik) / new_jj; *\/ */
                            /* out[j][k] = (out_jk - ratio*out_ik) * new_jj_inv; */
                            //out[j][j] = out[j][j] / new_mat[j][j];
			    out[j][k] = out_jk_val * new_jj_inv;
			    /* out[j][k] = out_jk_val / new_jj; */
                        }
                        else{

			    /* word_t new_jk = new_mat[j][k]; */
			    /* word_t new_ik = new_mat[i][k]; */
			    /* word_t new_jk_val = new_jk - ratio*new_ik; */
                            new_mat[j][k] = new_jk_val;

			    out[j][k] = out_jk_val;

                            //if(i > 0){
			    /* word_t out_jk = out[j][k]; */
			    /* word_t out_ik = out[i][k]; */
			    /* out[j][k] = out_jk - ratio*out_ik; */
		            /* } */
                            /* else{ //(i == 0) */

                            /*     if(i == k) */
			    /* 	{ */
                            /*         out[i][k] = 1; */
                            /*         out[j][k] = 0 - ratio; */
		            /*     } */
                            /*     else */
			    /* 	{ */
                            /*         out[i][k] = 0; */
			    /* 	    if(j == k){ */
			    /*             /\* if(i == k) *\/ */
			    /*             /\*     out[j][k] = 1 - ratio; *\/ */
                            /*             /\* else *\/ */
			    /*             out[j][k] = 1; */
                            /*             //out[j][k] = 1 - ratio*out[i][k]; */
		            /*         } */
			    /* 	    else */
			    /*             out[j][k] = 0; */

                                /* if(j == k){ */
                                /*     /\* if(i == k) *\/ */
                                /*     /\*     out[j][k] = 1 - ratio; *\/ */
                                /*     /\* else *\/ */
			        /*     out[j][k] = 1; */
                                /*     //out[j][k] = 1 - ratio*out[i][k]; */
                                /* } */
                                /* else{ */
                                /*     if(i == k) */
                                /*         out[j][k] = 0 - ratio; */
                                /*     else */
                                /*         out[j][k] = 0; */
                                /*     //out[j][k] = 0 - ratio*out[i][k]; */
                                /* } */
		        }
		    }
		}
            }
        }



	/* Row Operation to Make Principal Diagonal to 1 */
        loop_3:for(i = 0; i < DIM; i++)
        {
        #pragma HLS PIPELINE
	    word_t out_i = out[DIM-1][i];
	    word_t new_last = new_mat[DIM-1][DIM-1];
            out[DIM-1][i] = out_i / new_last;
        }

    }

    return 0;
}

/*
 * inverse.h
 *
 *  Created on: Aug 14, 2023
 *      Author: geichler
 */

#ifndef INVERSE_H_
#define INVERSE_H_


#include <stdio.h>
#include <stdlib.h>
//#include <conio.h>
#include <math.h>
#include <time.h>

#include "KalmanFilter.h"
#include "matmul.h"

//#define SIZE 1000

void inverse(float mat[N][N], float out[N][N])
{
		 //float a[SIZE][SIZE], x[SIZE], ratio;
		 float a[SIZE][SIZE], ratio;
		 int i,j,k,n;
		 time_t t;
		 //clrscr();

		 srand((unsigned) time(&t));

		 /* Inputs */
		 /* 1. Reading order of matrix */
		 printf("Enter order of matrix: ");
//		 scanf("%d", &n);

		 n = N;
		 printf("%d \n", n);

		 /* 2. Reading Matrix */
		 printf("Enter coefficients of Matrix:\n");
		 for(i=1;i<=n;i++)
		 {
			  for(j=1;j<=n;j++)
			  {
				   //printf("a[%d][%d] = ",i,j);
				   //scanf("%f", &a[i][j]);
				   a[i][j] = rand() % 50;
				   //printf("%f\n", a[i][j]);
				   mat[i-1][j-1] = (float)a[i][j];
			  }
		 }
		 /* Augmenting Identity Matrix of Order n */
		 for(i=1;i<=n;i++)
		 {
			  for(j=1;j<=n;j++)
			  {
				   if(i==j)
				   {
				    	a[i][j+n] = 1;
				   }
				   else
				   {
				    	a[i][j+n] = 0;
				   }
			  }
		 }
		 /* Applying Gauss Jordan Elimination */
		 for(i=1;i<=n;i++)
		 {
			  if(a[i][i] == 0.0)
			  {
				   printf("Mathematical Error!");
				   exit(0);
			  }
			  for(j=1;j<=n;j++)
			  {
				   if(i!=j)
				   {
					    ratio = a[j][i]/a[i][i];
					    for(k=1;k<=2*n;k++)
					    {
					     	a[j][k] = a[j][k] - ratio*a[i][k];
					    }
				   }
			  }
		 }
		 /* Row Operation to Make Principal Diagonal to 1 */
		 for(i=1;i<=n;i++)
		 {
			  for(j=n+1;j<=2*n;j++)
			  {
			   	a[i][j] = a[i][j]/a[i][i];
			  }
		 }
		 /* Displaying Inverse Matrix */
		 //printf("\nInverse Matrix is:\n");
		 for(i=1;i<=n;i++)
		 {
			  for(j=n+1;j<=2*n;j++)
			  {
			   	//printf("%f\t",a[i][j]);
			   	out[i-1][j-1-n] = a[i][j];
			  }
			  //printf("\n");
		 }
		 //getch();
		 return;
}

void inverse_clean(float new_mat[N][N], float out[N][N])
{

		 float ratio;
		 int i,j,k;

		 /* Applying Gauss Jordan Elimination */
		 for(i = 0; i < N; i++)
		 {
			  for(j = 0; j < N; j++)
			  {
				   if(i != j)
				   {
					    ratio = new_mat[j][i]/new_mat[i][i];
					    for(k = 0; k < N; k++)
					    {

					    	if(i == N-1){
					    		if(k == 0){//Calc the diagonal element first
					    			new_mat[j][j] = new_mat[j][j] - ratio*new_mat[i][j];
					    			//out[j][j] = (out[j][j] - ratio*out[i][j]) / new_mat[j][j];
					    		}
					    		else if(k == j){
					    			new_mat[j][0] = new_mat[j][0] - ratio*new_mat[i][0];
					    			//out[j][0] = (out[j][0] - ratio*out[i][0]) / new_mat[j][j];
					    		}
					    		else{
					    			new_mat[j][k] = new_mat[j][k] - ratio*new_mat[i][k];
					    			//out[j][k] = (out[j][k] - ratio*out[i][k]) / new_mat[j][j];
					    		}

					    		out[j][k] = (out[j][k] - ratio*out[i][k]) / new_mat[j][j];
					    		//out[j][j] = out[j][j] / new_mat[j][j];
					    	}
					    	else{

								new_mat[j][k] = new_mat[j][k] - ratio*new_mat[i][k];

								if(i > 0)
									out[j][k] = out[j][k] - ratio*out[i][k];
								else{ //(i == 0)
									if(i == k)
										out[i][k] = 1;
									else
										out[i][k] = 0;
									if(j == k){
										if(i == k)
											out[j][k] = 1 - ratio;
										else
											out[j][k] = 1;
										//out[j][k] = 1 - ratio*out[i][k];
									}
									else{
										if(i == k)
											out[j][k] = 0 - ratio;
										else
											out[j][k] = 0;
										//out[j][k] = 0 - ratio*out[i][k];
									}
								}
					    	}

					    }

				   }
//				   else{//(i == j)
//					   out[i][j] = out[i][j] / new_mat[i][i];
//				   }
			  }
		 }

		 /* Row Operation to Make Principal Diagonal to 1 */
		 for(i = 0; i < N; i++)
		 {
			 out[N-1][i] = out[N-1][i] / new_mat[N-1][N-1];
		 }

//		 /* Row Operation to Make Principal Diagonal to 1 */
//		 for(i = 0; i < N; i++)
//		 {
//			  for(j = 0; j < N;j++)
//			  {
//			   	out[i][j] = out[i][j]/new_mat[i][i];
//			  }
//		 }

		 return;
}

template<int DIM>
void inverse_clean(float new_mat[DIM][DIM], float out[DIM][DIM])
{

		 float ratio;
		 int i,j,k;

		 if(DIM == 2){
			 float a = new_mat[0][0];
			 float b = new_mat[0][1];
			 float c = new_mat[1][0];
			 float d = new_mat[1][1];

			 float det = (a * d) - (b * c);

			 if (det == 0) {
			     printf("The matrix is not invertible.\n");
			     return;
			 }

			 out[0][0] = d / det;
			 out[0][1] = (-1) * b / det;
			 out[1][0] = (-1) * c / det;
			 out[1][1] = a / det;

			 return;
		 }

		 /* Applying Gauss Jordan Elimination */
		 for(i = 0; i < DIM; i++)
		 {
			  for(j = 0; j < DIM; j++)
			  {
				   if(i != j)
				   {
					    ratio = new_mat[j][i]/new_mat[i][i];
					    for(k = 0; k < DIM; k++)
					    {

					    	if(i == DIM-1){
					    		if(k == 0){//Calc the diagonal element first
					    			new_mat[j][j] = new_mat[j][j] - ratio*new_mat[i][j];
					    			//out[j][j] = (out[j][j] - ratio*out[i][j]) / new_mat[j][j];
					    		}
					    		else if(k == j){
					    			new_mat[j][0] = new_mat[j][0] - ratio*new_mat[i][0];
					    			//out[j][0] = (out[j][0] - ratio*out[i][0]) / new_mat[j][j];
					    		}
					    		else{
					    			new_mat[j][k] = new_mat[j][k] - ratio*new_mat[i][k];
					    			//out[j][k] = (out[j][k] - ratio*out[i][k]) / new_mat[j][j];
					    		}

					    		out[j][k] = (out[j][k] - ratio*out[i][k]) / new_mat[j][j];
					    		//out[j][j] = out[j][j] / new_mat[j][j];
					    	}
					    	else{

								new_mat[j][k] = new_mat[j][k] - ratio*new_mat[i][k];

								if(i > 0)
									out[j][k] = out[j][k] - ratio*out[i][k];
								else{ //(i == 0)
									if(i == k)
										out[i][k] = 1;
									else
										out[i][k] = 0;
									if(j == k){
										if(i == k)
											out[j][k] = 1 - ratio;
										else
											out[j][k] = 1;
										//out[j][k] = 1 - ratio*out[i][k];
									}
									else{
										if(i == k)
											out[j][k] = 0 - ratio;
										else
											out[j][k] = 0;
										//out[j][k] = 0 - ratio*out[i][k];
									}
								}
					    	}

					    }

				   }
//				   else{//(i == j)
//					   out[i][j] = out[i][j] / new_mat[i][i];
//				   }
			  }
		 }

		 /* Row Operation to Make Principal Diagonal to 1 */
		 for(i = 0; i < DIM; i++)
		 {
			 out[DIM-1][i] = out[DIM-1][i] / new_mat[DIM-1][DIM-1];
		 }

//		 /* Row Operation to Make Principal Diagonal to 1 */
//		 for(i = 0; i < N; i++)
//		 {
//			  for(j = 0; j < N;j++)
//			  {
//			   	out[i][j] = out[i][j]/new_mat[i][i];
//			  }
//		 }

		 return;
}

template<int DIM>
void iterative_inverse(float new_mat[DIM][DIM], float out[DIM][DIM], float new_out_final[DIM][DIM])
{
	float new_out[DIM][DIM] = {{0}};
	//float new_out_final[DIM][DIM] = {{0}};
	//float I2[DIM][DIM];

	//A*X_n
	matmul<DIM>(new_mat, out, new_out);

	//2I - A*X_n
	for(int i = 0; i < DIM; i++)
		for(int j = 0; j < DIM; j++)
			if(i == j)
				new_out[i][j] = 2 - new_out[i][j];
			else
				new_out[i][j] = 0 - new_out[i][j];

	//X_n*(2I-A*X_n)
	matmul<DIM>(out, new_out, new_out_final);
}

template<int DIM>
void many_iterative_inverse(float new_mat[DIM][DIM], float out[DIM][DIM], float new_out_final[DIM][DIM], int iter)
{
	for(int k = 0; k < iter; k++){

		//std::cout << "ITERATION " << k << std::endl;

		//float next_out[DIM][DIM] = {{0}};
		iterative_inverse<DIM>(new_mat, out, new_out_final);
		//iterative_inverse<N>(mat, rand_inv, next_out);

		//Initialize inverse as the new iteration
		for(int i = 0; i < DIM; i++)
			for(int j = 0; j < DIM; j++)
				out[i][j] = new_out_final[i][j];
				//rand_inv[i][j] = next_out[i][j];

	}
}

template<int DIM>
void inverse_free(float new_mat[DIM][DIM], float out[DIM][DIM])
{

	for(int i = 0; i < DIM; i++)
	{
		for(int j = 0; j < DIM; j++)
		{
			if(i == j) //Diagonal
			{
				out[i][j] = 1 / new_mat[i][j];
			}
			else
			{
				out[i][j] = /*pow(-1, i + j)*/ (-1) * (1 / new_mat[i][i]) * new_mat[i][j] * (1 / new_mat[j][j]);
			}
		}
	}
}

//template<int X_DIM, int Z_DIM>
//void compute_S_inverse_free(float S_array[Z_DIM][Z_DIM], float P_array[X_DIM][X_DIM], float H_array[Z_DIM][X_DIM], float R_array[Z_DIM][Z_DIM])
//{
//
//	for(int i = 0; i < Z_DIM; i++)
//	{
//		for(int j = 0; j < Z_DIM; j++)
//		{
//			if(i == j) //Diagonal
//			{
//				out[i][j] = 1 / new_mat[i][j];
//			}
//			else
//			{
//				out[i][j] = pow(-1, i + j) * (1 / new_mat[i][i]) * new_mat[i][j] * (1 / new_mat[j][j]);
//			}
//		}
//	}
//}

#endif /* INVERSE_H_ */

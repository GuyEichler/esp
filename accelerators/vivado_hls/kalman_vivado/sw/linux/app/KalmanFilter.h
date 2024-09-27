/*
 * KalmanFilter.h
 *
 *  Created on: Aug 10, 2023
 *      Author: geichler
 */

#ifndef KALMANFILTER_H_
#define KALMANFILTER_H_

#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <bits/stdc++.h>
#include <time.h>
//#include <omp.h>

#include <eigen/Eigen/Dense>
#include <eigen/Eigen/SVD>
#include <algorithm>
#include <vector>

#define N 10
#define TINY 1.0e-20
#define SIZE 1000


#include "A_array.h"
#include "H_array.h"
#include "W_array.h"
#include "Q_array.h"
#include "initial_state_array.h"
#include "measurements_array.h"
#include "prediction_array.h"
#include "real_array.h"

////#include "A_array_soma.h"
////#include "H_array_soma.h"
////#include "W_array_soma.h"
////#include "Q_array_soma.h"
////#include "initial_state_array_soma.h"
////#include "measurements_array_soma.h"
////#include "prediction_array_soma.h"
////#include "real_array_soma.h"
//
////#include "A_array_hc.h"
////#include "H_array_hc.h"
////#include "W_array_hc.h"
////#include "Q_array_hc.h"
////#include "initial_state_array_hc.h"
////#include "measurements_array_hc.h"
////#include "prediction_array_hc.h"
////#include "real_array_hc.h"


#include "inverse.h"


using namespace Eigen;

namespace kf
{
	template<int DIM_X/*state vector*/, int DIM_Z/*measurement vectoe*/>
	class KalmanFilter
	{
	public:
		KalmanFilter() = default;

		void prediction(const Matrix<float, DIM_X, DIM_X> &Mat_F/*state transition matrix*/,
				const Matrix<float, DIM_X, DIM_X> &Mat_Q/*process noise covariance*/)
		{
			mVec_X = Mat_F * mVec_X;
			mMat_P = Mat_F * mMat_P * Mat_F.transpose() + Mat_Q;
		}

		void correction(const Vector<float, DIM_Z> &Vec_Z/*measurement vector*/,
				const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/)
		{
			//Innovation / residual
			Vector<float, DIM_Z> Vec_Y = Vec_Z - (Mat_H * mVec_X);

			//Innovation covariance
			Matrix<float, DIM_Z, DIM_Z> Mat_S = Mat_H * mMat_P * Mat_H.transpose() + Mat_R;

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K = mMat_P * Mat_H.transpose() * Mat_S.inverse();

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			//Update state vector
			mVec_X = mVec_X  + Mat_K * Vec_Y;

			//Update state covariance
			mMat_P = (Mat_I -  Mat_K * Mat_H) * mMat_P;

		}

		Vector<float, DIM_X>& Vec_X() {return mVec_X;}
		const Vector<float, DIM_X>& Vec_X() const {return mVec_X;}

		Matrix<float, DIM_X, DIM_X>& Mat_P() {return mMat_P;}
		const Matrix<float, DIM_X, DIM_X>& Mat_P() const {return mMat_P;}

	private:
		//State vector
		Vector<float, DIM_X> mVec_X;

		//State covariance matrix
		Matrix<float, DIM_X, DIM_X> mMat_P;

	};

	template<int DIM_X/*state vector*/, int DIM_Z/*measurement vector*/>
	class KalmanFilterRegression
	{
	public:
		KalmanFilterRegression() = default;

		void fit(float scale//const MatrixXf &Mat_X, const MatrixXf &Mat_Y
				)
		{
			MatrixXf X = mMat_Ytrain.transpose();
			MatrixXf Z = mMat_Xtrain.transpose();

			float nt = (float)X.cols();

			MatrixXf X2 = X.block(0,1,X.rows(),X.cols()-1);
			MatrixXf X1 = X.block(0,0,X.rows(),X.cols()-2);

			MatrixXf X11 = X1 * X1.transpose();

			mMat_F = X2 * X1.transpose() * X11.inverse();

			mMat_Q = (X2 - mMat_F * X1) * (X2 - mMat_F * X1).transpose()/(nt-1)/scale;

			mMat_H = Z * X.transpose() * (X * X.transpose()).inverse();

			mMat_R = ((Z - mMat_H * X)*((Z - mMat_H * X).trnaspose())) / nt;
		}

		void prediction(//const Matrix<float, DIM_X, DIM_X> &Mat_F/*state transition matrix*/,
				//const Matrix<float, DIM_X, DIM_X> &Mat_Q/*process noise covariance*/
				)
		{
			mVec_X = mMat_F * mVec_X;
			mMat_P = mMat_F * mMat_P * mMat_F.transpose() + mMat_Q;
		}

		void correction(const Vector<float, DIM_Z> &Vec_Z/*measurement vector*/
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Innovation / residual
			Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mVec_X);

			//Innovation covariance
			Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * mMat_P * mMat_H.transpose() + mMat_R;

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K = mMat_P * mMat_H.transpose() * Mat_S.inverse();

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			//Update state vector
			mVec_X = mVec_X  + Mat_K * Vec_Y;

			//Update state covariance
			mMat_P = (Mat_I -  Mat_K * mMat_H) * mMat_P;

		}

		void predict_and_correct(const Vector<float, DIM_Z> &Vec_Z/*measurement vector*/
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Merge prediction into the computations of correction

			//Innovation / residual
			Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mMat_F * mVec_X);

			//Innovation covariance
			Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() + mMat_R;

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K = (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() * Mat_S.inverse();

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			//Update state vector
			mVec_X = mMat_F * mVec_X  + Mat_K * Vec_Y;

			//Update state covariance
			mMat_P = (Mat_I -  Mat_K * mMat_H) * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q);
			//std::cout << "Matrix P = " << mMat_P << std::endl;

		}

		void predict_and_correct_opt(const Vector<float, DIM_Z> &Vec_Z/*measurement vector*/
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Merge prediction into the computations of correction

			//Innovation / residual
			//Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mMat_F * mVec_X);

			//Innovation covariance
			//Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() + mMat_R;

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K = (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() * (mMat_H * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() + mMat_R).inverse();

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			//Update state vector
			mVec_X = mMat_F * mVec_X  + Mat_K * (Vec_Z - (mMat_H * mMat_F * mVec_X));

			//Update state covariance
			mMat_P = (Mat_I -  Mat_K * mMat_H) * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q);

		}

		void correction_inv(const Vector<float, DIM_Z> &Vec_Z/*measurement vector*/
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Innovation / residual
			Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mVec_X);

			//Innovation covariance
			Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * mMat_P * mMat_H.transpose() + mMat_R;

			//Invert the S matrix manually
			Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_row = Mat_S;
			//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_inv;
			float S_array[DIM_Z][DIM_Z];
			float inv_array[DIM_Z][DIM_Z];
			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row;
			inverse_clean<DIM_Z>(S_array, inv_array);

			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&inv_array[0][0], DIM_Z, DIM_Z);

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K = mMat_P * mMat_H.transpose() * Mat_S_inv;

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			//Update state vector
			mVec_X = mVec_X  + Mat_K * Vec_Y;

			//Update state covariance
			mMat_P = (Mat_I -  Mat_K * mMat_H) * mMat_P;

		}

		void correction_inv_iter(const Vector<float, DIM_Z> &Vec_Z/*measurement vector*/, int iter
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Innovation / residual
			Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mVec_X);

			//Innovation covariance
			Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * mMat_P * mMat_H.transpose() + mMat_R;

			//Invert the S matrix manually
			Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_row = Mat_S;
			//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_inv;
			float S_array[DIM_Z][DIM_Z];
			float inv_array[DIM_Z][DIM_Z], final_inv_array[DIM_Z][DIM_Z];
			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row;
			inverse_clean<DIM_Z>(S_array, inv_array);

			many_iterative_inverse<DIM_Z>(S_array, inv_array, final_inv_array, iter);

			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&final_inv_array[0][0], DIM_Z, DIM_Z);

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K = mMat_P * mMat_H.transpose() * Mat_S_inv;

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			//Update state vector
			mVec_X = mVec_X  + Mat_K * Vec_Y;

			//Update state covariance
			mMat_P = (Mat_I -  Mat_K * mMat_H) * mMat_P;

		}

		void predict_and_correct_inv(const Vector<float, DIM_Z> &Vec_Z/*measurement vector*/
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Merge prediction into the computations of correction

			//Innovation / residual
			Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mMat_F * mVec_X);

			//Innovation covariance
			Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() + mMat_R;


			//Invert the S matrix manually
			Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_row = Mat_S;
			//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_inv;
			float S_array[DIM_Z][DIM_Z];
			float (*inv_array)[DIM_Z];//[DIM_Z];
			inv_array = mInv_array1;

			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row;
			inverse_clean<DIM_Z>(S_array, inv_array);

			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&inv_array[0][0], DIM_Z, DIM_Z);

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K = (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() * Mat_S_inv;

			//std::cout << Mat_K << std::endl;

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			//Update state vector
			mVec_X = mMat_F * mVec_X  + Mat_K * Vec_Y;

			//Update state covariance
			mMat_P = (Mat_I -  Mat_K * mMat_H) * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q);

		}

		void predict_and_correct_inv_manual(//const Vector<float, DIM_Z> &Vec_Z_eigen,/*measurement vector*/
				float Vec_Z[DIM_Z][1]
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Merge prediction into the computations of correction

//			Map<Matrix<float, DIM_X, 1>> printX(&Vec_X_arr[0][0]);
//			std::cout << printX << std::endl;

			//Innovation / residual
//			Vector<float, DIM_Z> Vec_Y_eigen = Vec_Z_eigen - (mMat_H * mMat_F * mVec_X);
			float Vec_Y[DIM_Z][1];
			float tmp_zx[DIM_Z][DIM_X];
			matmul<DIM_Z, DIM_X, DIM_X>(Mat_H_arr, Mat_F_arr, tmp_zx);
			matmul<DIM_Z, DIM_X, 1>(tmp_zx, Vec_X_arr, Vec_Y);
			matsub<DIM_Z, 1>(Vec_Z, Vec_Y, Vec_Y);

//			Map<Matrix<float, DIM_Z, 1>> printZ(&Vec_Z[0][0]);
//			std::cout << printZ << std::endl;
//			std::cout << std::endl;
//			std::cout << std::endl;
//			std::cout << Vec_Z_eigen << std::endl;

			//Innovation covariance
//			Matrix<float, DIM_Z, DIM_Z> Mat_S_eigen = mMat_H * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() + mMat_R;
			float Mat_S[DIM_Z][DIM_Z];
			float tmp_xx1[DIM_X][DIM_X];
			float tmp_xx2[DIM_X][DIM_X];
			float tmp_xz[DIM_X][DIM_Z];
			matmul<DIM_X, DIM_X, DIM_X>(Mat_F_arr, Mat_P_arr, tmp_xx1);
			matmul_trans<DIM_X, DIM_X, DIM_X>(tmp_xx1, Mat_F_arr, tmp_xx2);
			matadd<DIM_X, DIM_X>(tmp_xx2, Mat_Q_arr, tmp_xx1);
			matmul_trans<DIM_X, DIM_X, DIM_Z>(tmp_xx1, Mat_H_arr, tmp_xz);
			matmul<DIM_Z, DIM_X, DIM_Z>(Mat_H_arr, tmp_xz, Mat_S);
			matadd<DIM_Z, DIM_Z>(Mat_S, Mat_R_arr, Mat_S);

//			Map<Matrix<float, DIM_Z, DIM_Z>> printS(&Mat_S[0][0]);
//			std::cout << printS << std::endl;

			//Invert the S matrix manually
//			Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_row = Mat_S_eigen;
//			//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_inv;
//			float S_array[DIM_Z][DIM_Z];
			float (*S_inv)[DIM_Z];//[DIM_Z];
			S_inv = mInv_array1;

//			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row;
//			inverse_clean<DIM_Z>(S_array, S_inv);
//			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&S_inv[0][0], DIM_Z, DIM_Z);

			inverse_clean<DIM_Z>(Mat_S, S_inv);

//			Map<Matrix<float, DIM_Z, DIM_Z>> printSinv(&S_inv[0][0]);
//			std::cout << printSinv << std::endl;



			//Kalman gain
//			Matrix<float, DIM_X, DIM_Z> Mat_K_eigen = (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() * Mat_S_inv;
			float Mat_K[DIM_X][DIM_Z];
			matmul<DIM_X, DIM_Z, DIM_Z>(tmp_xz, S_inv, Mat_K);

//			Map<Matrix<float, DIM_X, DIM_Z, RowMajor>> printK(&Mat_K[0][0]);
//			std::cout << printK << std::endl;

			//Identity matrix
//			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();
			float Mat_I[DIM_X][DIM_X];
			for(int i = 0; i < DIM_X; i++)
				for(int j = 0; j < DIM_X; j++)
				{
					if(i == j) Mat_I[i][j] = 1.0;
					else Mat_I[i][j] = 0.0;
				}

			//Update state vector
//			Map<Matrix<float, DIM_Z, 1>> Vec_Y_eigen(&Vec_Y[0][0]);
//			Matrix<float, DIM_X, DIM_Z, RowMajor> Mat_K_eigen(&Mat_K[0][0]);
//			mVec_X = mMat_F * mVec_X  + Mat_K_eigen * Vec_Y_eigen;
			float tmp_x1[DIM_X][1];
			float tmp_x2[DIM_X][1];
			matmul<DIM_X, DIM_X, 1>(Mat_F_arr, Vec_X_arr, tmp_x1);
			matmul<DIM_X, DIM_Z, 1>(Mat_K, Vec_Y, tmp_x2);
			matadd<DIM_X, 1>(tmp_x2, tmp_x1, Vec_X_arr);

			/* Map<Matrix<float, DIM_X, 1>> Vec_X_eigen(&Vec_X_arr[0][0]); */
			/* mVec_X = Vec_X_eigen; */

//			Map<Matrix<float, DIM_X, 1>> printX1(&Vec_X_arr[0][0]);
//			std::cout << printX1 << std::endl;

			//Update state covariance
//			mMat_P = (Mat_I -  Mat_K_eigen * mMat_H) * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q);
			matmul<DIM_X, DIM_Z, DIM_X>(Mat_K, Mat_H_arr ,tmp_xx2);
			matsub<DIM_X, DIM_X>(Mat_I, tmp_xx2 ,tmp_xx2);
			matmul<DIM_X, DIM_X, DIM_X>(tmp_xx2, tmp_xx1, Mat_P_arr);


			/* Map<Matrix<float, DIM_X, DIM_X, RowMajor>> Mat_P_eigen(&Mat_P_arr[0][0]); */
			/* mMat_P = Mat_P_eigen; */

		}

		void predict_and_correct_iter(const Vector<float, DIM_Z> &Vec_Z, int iter/*measurement vector*/
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Merge prediction into the computations of correction

			//Innovation / residual
			Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mMat_F * mVec_X);

			//Innovation covariance
			Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() + mMat_R;

			//Invert the S matrix manually
			Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_row = Mat_S;
			//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_inv;
			float S_array[DIM_Z][DIM_Z];
			//float inv_array[DIM_Z][DIM_Z];
			float (*final_inv_array)[DIM_Z], (*inv_array)[DIM_Z]; //[DIM_Z][DIM_Z];

			if(mPingpong){
				inv_array = mInv_array1;
				final_inv_array = mInv_array2;
			}else{
				inv_array = mInv_array2;
				final_inv_array = mInv_array1;
			}

			mPingpong = !mPingpong;

			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row;
			//inverse_clean<DIM_Z>(S_array, inv_array);

			many_iterative_inverse<DIM_Z>(S_array, inv_array, final_inv_array, iter);

			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&final_inv_array[0][0], DIM_Z, DIM_Z);

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K = (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() * Mat_S_inv;

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			//Update state vector
			mVec_X = mMat_F * mVec_X  + Mat_K * Vec_Y;

			//Update state covariance
			mMat_P = (Mat_I -  Mat_K * mMat_H) * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q);

		}

		void predict_and_correct_sskf(const Vector<float, DIM_Z> &Vec_Z/*measurement vector*/
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Merge prediction into the computations of correction

			//Innovation / residual
			Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mMat_F * mVec_X);

			//Innovation covariance
			//Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() + mMat_R;

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K = mMat_SSK;

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			//Update state vector
			mVec_X = mMat_F * mVec_X  + Mat_K * Vec_Y;
			std::cout << "Vec X = " << mVec_X << std::endl;

			//Update state covariance
			mMat_P = (Mat_I -  mMat_SSK * mMat_H) * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q);
			std::cout << "Matrix P = " << mMat_P << std::endl;

		}

		void predict_and_correct_sskf_iter(const Vector<float, DIM_Z> &Vec_Z, int iter/*measurement vector*/
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Merge prediction into the computations of correction

			//Innovation / residual
			Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mMat_F * mVec_X);

			//Innovation covariance
			Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() + mMat_R;

			//Kalman gain
			//Matrix<float, DIM_X, DIM_Z> Mat_K = mMat_SSK;

			//Invert the S matrix manually
			//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_row = Mat_S;
			//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_inv;
			float S_array[DIM_Z][DIM_Z];
			float inv_array[DIM_Z][DIM_Z], final_inv_array[DIM_Z][DIM_Z];
			//Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row;
			//inverse_free<DIM_Z>(S_array, inv_array);

			for(int i = 0; i < DIM_Z; i++)
				for(int j = 0; j < DIM_Z; j++)
				{
					S_array[i][j] = Mat_S(i, j);
					inv_array[i][j] = mInv_array1[i][j];
				}

			//Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row;
			//inverse_clean<DIM_Z>(S_array, inv_array);

			many_iterative_inverse<DIM_Z>(S_array, inv_array, final_inv_array, iter);

			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&inv_array[0][0], DIM_Z, DIM_Z);
			//Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&inv_array[0][0], DIM_Z, DIM_Z);

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K = (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() * Mat_S_inv;

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			//Update state vector
			mVec_X = mMat_F * mVec_X  + Mat_K * Vec_Y;

			//Update state covariance
			mMat_P = (Mat_I -  mMat_SSK * mMat_H) * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q);
			//std::cout << "Matrix P = " << mMat_P << std::endl;

		}

		void compute_ssk_gain(int iterations//const Vector<float, DIM_Z> &Vec_Z/*measurement vector*/
						//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
						//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
						)
		{
			//Set mMat_SSK to the steady-state value
			Matrix<float, DIM_X, DIM_X> Mat_P = Matrix<float, DIM_X, DIM_X>::Zero();
			Matrix<float, DIM_Z, DIM_Z> Mat_S;
			for(int i = 0; i < iterations; i++)
			{
				Mat_S = mMat_H * (mMat_F * Mat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() + mMat_R;
				mMat_SSK = (mMat_F * Mat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() * Mat_S.inverse();
				Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();
				Mat_P = (Mat_I -  mMat_SSK * mMat_H) * (mMat_F * Mat_P * mMat_F.transpose() + mMat_Q);
			}

			//std::cout << "Matrix K = " << mMat_SSK << std::endl;

			//Store the converged inverse
			//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_inv_row = Mat_S.inverse();
			Matrix<float, DIM_Z, DIM_Z> Mat_S_inv = Mat_S.inverse();
			//float (*inv_array)[DIM_Z];//[DIM_Z];
			//inv_array = mInv_array1;
			//Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&mInv_array1[0][0], DIM_Z, DIM_Z) = Mat_S_inv_row;
			for(int i = 0; i < DIM_Z; i++)
				for(int j = 0; j < DIM_Z; j++)
						mInv_array1[i][j] = Mat_S_inv(i, j);
		}

		void predict_and_correct_ifkf(const Vector<float, DIM_Z> &Vec_Z/*measurement vector*/
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Merge prediction into the computations of correction

			//Innovation / residual
			Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mMat_F * mVec_X);

			Matrix<float, DIM_X, DIM_X> P_pred = (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q);
			//Matrix<float, DIM_X, DIM_X, RowMajor> Mat_P_pred_row = P_pred;

//			float P_array[DIM_X][DIM_X];
//			Map<Matrix<float, DIM_X, DIM_X, RowMajor>>(&P_array[0][0], DIM_Z, DIM_Z) = Mat_P_pred_row;
//
//			Matrix<float, DIM_Z, DIM_X, RowMajor> Mat_H_pred_row = mMat_H;
//			float H_array[DIM_Z][DIM_X];
//			Map<Matrix<float, DIM_Z, DIM_X, RowMajor>>(&H_array[0][0], DIM_Z, DIM_X) = Mat_H_pred_row;
//
//			Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_R_pred_row = mMat_R;
//			float R_array[DIM_Z][DIM_Z];
//			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&H_array[0][0], DIM_Z, DIM_Z) = Mat_R_pred_row;

			//Innovation covariance
			Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * P_pred * mMat_H.transpose() + mMat_R;

//			float S_array[DIM_Z][DIM_Z];

//			compute_S_inverse_free<DIM_X,DIM_Z>(S_array, P_array, H_array, R_array);

			//std::cout << mMat_R << std::endl;

//			std::cout << "REAL INVERSE" << std::endl;
//			std::cout << std::endl;
//			std::cout << std::endl;
//			std::cout << std::endl;
//			std::cout << Mat_S.inverse() << std::endl;
//			std::cout << std::endl;
//			std::cout << std::endl;
//			std::cout << std::endl;

			//Invert the S matrix manually
			Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_row = Mat_S;
			//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_inv;
			float S_array[DIM_Z][DIM_Z];
			float inv_array[DIM_Z][DIM_Z];
			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row;
			inverse_free<DIM_Z>(S_array, inv_array);

			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&inv_array[0][0], DIM_Z, DIM_Z);

//			std::cout << " ====================== IFKF INVERSE =================" << std::endl;
//			std::cout << std::endl;
//			std::cout << std::endl;
//			std::cout << std::endl;
//			std::cout << Mat_S_inv << std::endl;
//			std::cout << std::endl;
//			std::cout << std::endl;
//			std::cout << std::endl;

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K = P_pred * mMat_H.transpose() * Mat_S_inv;

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			//Update state vector
			mVec_X = mMat_F * mVec_X  + Mat_K * Vec_Y;

			//Update state covariance
			mMat_P = (Mat_I -  Mat_K * mMat_H) * P_pred;

		}

		void predict_and_correct_ifkf_iter(const Vector<float, DIM_Z> &Vec_Z, int iter/*measurement vector*/
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Merge prediction into the computations of correction

			//Innovation / residual
			Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mMat_F * mVec_X);

			//Innovation covariance
			Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() + mMat_R;

			//Invert the S matrix manually
			Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_row = Mat_S;
			//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_inv;
			float S_array[DIM_Z][DIM_Z];
			float inv_array[DIM_Z][DIM_Z], final_inv_array[DIM_Z][DIM_Z];
			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row;
			inverse_free<DIM_Z>(S_array, inv_array);

			//Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row;
			//inverse_clean<DIM_Z>(S_array, inv_array);

			many_iterative_inverse<DIM_Z>(S_array, inv_array, final_inv_array, iter);

			//Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&final_inv_array[0][0], DIM_Z, DIM_Z);
			Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&inv_array[0][0], DIM_Z, DIM_Z);

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K = (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() * Mat_S_inv;

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			//Update state vector
			mVec_X = mMat_F * mVec_X  + Mat_K * Vec_Y;

			//Update state covariance
			mMat_P = (Mat_I -  Mat_K * mMat_H) * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q);

		}

		void predict_and_correct_taylor(const Vector<float, DIM_Z> &Vec_Z/*measurement vector*/
				//const Matrix<float, DIM_Z, DIM_X> &Mat_H/*measurement model*/,
				//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/*measurement noise covariance*/
				)
		{
			//Merge prediction into the computations of correction

			Matrix<float, DIM_X, DIM_X> Mat_P_pred = mMat_F * mMat_P * mMat_F.transpose() + mMat_Q;

			//Identity matrix
			Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity();

			Matrix<float, DIM_X, DIM_X> Mat_P_I = Mat_P_pred - Mat_I;

			//Calculate appoximate Kalman gain
			//compute_C();

			Matrix<float, DIM_Z, DIM_Z> Sum_Mat_term = Matrix<float, DIM_Z, DIM_Z>::Zero();

			for(int u = 0; u < DIM_X; u++)
				for(int l = 0; l < DIM_X; l++)
				{
					Matrix<float, DIM_Z, DIM_Z> Mat_term = compute_taylor_term(u, l, Mat_P_I(u, l));
					Sum_Mat_term += Mat_term;
				}

			Sum_Mat_term = mMat_C - Sum_Mat_term;

			//std::cout << "Sum = " <<  Sum_Mat_term << std::endl;

			//Kalman gain
			Matrix<float, DIM_X, DIM_Z> Mat_K_approx = Mat_P_pred * mMat_H.transpose() * Sum_Mat_term;

			//Innovation / residual
			Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mMat_F * mVec_X);

			//Innovation covariance
			//Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() + mMat_R;


			//Update state vector
			mVec_X = mMat_F * mVec_X  + Mat_K_approx * Vec_Y;

			//std::cout << "Vec X = " << mVec_X << std::endl;

			//Update state covariance
			mMat_P = (Mat_I -  Mat_K_approx * mMat_H) * Mat_P_pred * (Mat_I - Mat_K_approx * mMat_H).transpose()
						+ Mat_K_approx * mMat_R * Mat_K_approx.transpose();
			//std::cout << "Matrix P = " << mMat_P << std::endl;

		}

		void compute_C()
		{
			mMat_C = (mMat_H * mMat_H.transpose() + mMat_R).inverse();
			//std::cout << "Mat C = " << mMat_C << std::endl;
		}

		Matrix<float, DIM_Z, DIM_Z> compute_taylor_term(int u, int l, float scalar)
		{
			Matrix<float, DIM_X, DIM_X> Mat_E = Matrix<float, DIM_X, DIM_X>::Zero();
			Mat_E(u, l) = 1.0;

			Matrix<float, DIM_Z, DIM_Z> Mat_term = mMat_C * mMat_H * Mat_E * mMat_H.transpose() * mMat_C;

			//std::cout << "Mat_term = " << mMat_C * mMat_H * Mat_E << std::endl;

			//std::cout << "scalar (" << u << "," << l << ") = " << scalar << std::endl;

			return Mat_term * scalar;
		}

		void reset()
		{
			mVec_X = Vector<float, DIM_X>::Zero();
			mMat_P = Matrix<float, DIM_X, DIM_X>::Zero();

			//State transition matrix - A in neural decoding
			mMat_F = Matrix<float, DIM_X, DIM_X>::Zero();

			//Process noise covariance matrix - W in neural decoding
			mMat_Q = Matrix<float, DIM_X, DIM_X>::Zero();

			//Measurement noise covariance matrix - Q in neural decoding
			mMat_R = Matrix<float, DIM_Z, DIM_Z>::Zero();

			//Measurement model matrix - H in neural decoding
			mMat_H = Matrix<float, DIM_Z, DIM_X>::Zero();

			//Steady-state Kalman gain
			mMat_SSK = Matrix<float, DIM_X, DIM_Z>::Zero();

			//Taylor approximation Kalman gain
			mMat_taylor_K = Matrix<float, DIM_X, DIM_Z>::Zero();

			mPingpong = true;
		}

		void print_SSK_inv()
		{
			printf("Starting to write SSK inverse into a file...");

			std::ofstream myfile ("SSK_inv.h");
			if (myfile.is_open())
			{
				myfile << "float SSK_inv[] = {";
			}

			//int iter = time_stamps;
			for(int i = 0; i < DIM_Z; i++)
			{
				for(int j = 0; j < DIM_Z; j++)
				{
//				Matrix<float, neurons, neurons> mat_S = mat_H * (mat_F * mat_P * mat_F.transpose() + mat_Q) * mat_H.transpose() + mat_R;
//				Matrix<float, states, neurons> mat_K = (mat_F * mat_P * mat_F.transpose() + mat_Q) * mat_H.transpose() * mat_S.inverse();
//				mat_P = (mat_I -  mat_K * mat_H) * (mat_F * mat_P * mat_F.transpose() + mat_Q);
//
//				float P_flat[states*states];
//				Map<Matrix<float, states, states, RowMajor> >(&P_flat[0], mat_P.rows(), mat_P.cols()) = mat_P;

					if (myfile.is_open())
					{
						if(i == DIM_Z-1 && j == DIM_Z - 1)
							myfile << mInv_array1[i][j] << "}; ";
						else
							myfile << mInv_array1[i][j] << ", ";
					}
				}
			}

			myfile.close();

			printf("Finished writing SSK inverse into a file...\n");
		}

		void print_SSK_gain()
		{
			printf("Starting to write SSK gain into a file...");

			std::ofstream myfile ("SSK_gain.h");
			if (myfile.is_open())
			{
				myfile << "float SSK_gain[] = {";
			}

			//int iter = time_stamps;
			for(int i = 0; i < DIM_X; i++)
			{
				for(int j = 0; j < DIM_Z; j++)
				{
//				Matrix<float, neurons, neurons> mat_S = mat_H * (mat_F * mat_P * mat_F.transpose() + mat_Q) * mat_H.transpose() + mat_R;
//				Matrix<float, states, neurons> mat_K = (mat_F * mat_P * mat_F.transpose() + mat_Q) * mat_H.transpose() * mat_S.inverse();
//				mat_P = (mat_I -  mat_K * mat_H) * (mat_F * mat_P * mat_F.transpose() + mat_Q);
//
//				float P_flat[states*states];
//				Map<Matrix<float, states, states, RowMajor> >(&P_flat[0], mat_P.rows(), mat_P.cols()) = mat_P;

					if (myfile.is_open())
					{
						if(i == DIM_X-1 && j == DIM_Z - 1)
							myfile << mMat_SSK(i,j) << "}; ";
						else
							myfile << mMat_SSK(i,j) << ", ";
					}
				}
			}

			myfile.close();

			printf("Finished writing SSK inverse into a file...\n");
		}

		void print_C()
		{
			printf("Starting to write C into a file...");

			std::ofstream myfile ("C.h");
			if (myfile.is_open())
			{
				myfile << "float C[] = {";
			}

			//int iter = time_stamps;
			for(int i = 0; i < DIM_Z; i++)
			{
				for(int j = 0; j < DIM_Z; j++)
				{
//				Matrix<float, neurons, neurons> mat_S = mat_H * (mat_F * mat_P * mat_F.transpose() + mat_Q) * mat_H.transpose() + mat_R;
//				Matrix<float, states, neurons> mat_K = (mat_F * mat_P * mat_F.transpose() + mat_Q) * mat_H.transpose() * mat_S.inverse();
//				mat_P = (mat_I -  mat_K * mat_H) * (mat_F * mat_P * mat_F.transpose() + mat_Q);
//
//				float P_flat[states*states];
//				Map<Matrix<float, states, states, RowMajor> >(&P_flat[0], mat_P.rows(), mat_P.cols()) = mat_P;

					if (myfile.is_open())
					{
						if(i == DIM_Z-1 && j == DIM_Z - 1)
							myfile << mMat_C(i,j) << "}; ";
						else
							myfile << mMat_C(i,j) << ", ";
					}
				}
			}

			myfile.close();

			printf("Finished writing C into a file...\n");
		}

		// Function to compute the cross-correlation matrix between measurements and states
		Matrix<float, DIM_Z, DIM_X> computeCrossCorrelationMatrix(const MatrixXf& measurements, const MatrixXf& states)
		{
			int numSamples = measurements.rows();

		    // Center the data (subtract mean)
		    MatrixXf centeredMeasurements = measurements.rowwise() - measurements.colwise().mean();
		    MatrixXf centeredStates = states.rowwise() - states.colwise().mean();

//		    std::cout << "Num of rows centeredMeasurements " << centeredMeasurements.rows() << std::endl;
//		    std::cout << "Num of cols centeredMeasurements " << centeredMeasurements.cols() << std::endl;
//		    std::cout << "Num of rows centeredStates " << centeredStates.rows() << std::endl;
//		    std::cout << "Num of cols centeredStates " << centeredStates.cols() << std::endl;

		    // Compute the covariance matrix between measurements and states
		    MatrixXf crossCovariance = (centeredMeasurements.transpose() * centeredStates) / float(numSamples - 1);

//		    std::cout << "Num of rows crossCovariance " << crossCovariance.rows() << std::endl;
//		    std::cout << "Num of cols crossCovariance " << crossCovariance.cols() << std::endl;

		    // Compute the standard deviations of measurements and states
		    VectorXf stdMeasurements = (centeredMeasurements.array().square().colwise().sum() / float(numSamples - 1)).sqrt();
		    VectorXf stdStates = (centeredStates.array().square().colwise().sum() / float(numSamples - 1)).sqrt();

//		    std::cout << "Size of stdMeasurements " << stdMeasurements.rows() << std::endl;
//		    std::cout << "Size of stdStates " << stdStates.rows() << std::endl;

		    // Normalize the cross-covariance matrix to get the cross-correlation matrix
		    Matrix<float, DIM_Z, DIM_X> crossCorrelation = crossCovariance.array().rowwise() / stdStates.transpose().array();
		    crossCorrelation = crossCorrelation.array().colwise() / stdMeasurements.array();

//		    std::cout << "Num of rows crossCorrelation " << crossCorrelation.rows() << std::endl;

		    return crossCorrelation;
		}

		std::vector<int> selectTopCorrelatedMeasurements(const Matrix<float, DIM_Z, DIM_X>& crossCorrelation, int numToSelect)
		{
			std::vector<int> selectedIndices;
		    int numMeasurements = crossCorrelation.rows();

		    // Compute the maximum correlation for each measurement
		    VectorXf maxCorrelation = crossCorrelation.cwiseAbs().rowwise().maxCoeff();

		    // Create a vector of indices
		    std::vector<int> indices(numMeasurements);
		    std::iota(indices.begin(), indices.end(), 0);  // Fill with 0, 1, 2, ..., numMeasurements-1

		    // Sort indices based on maximum correlation values in descending order
		    std::sort(indices.begin(), indices.end(), [&maxCorrelation](int i, int j) {
		        return maxCorrelation(i) > maxCorrelation(j);
		    });

		    // Select the top correlated measurements
		    selectedIndices.assign(indices.begin(), indices.begin() + numToSelect);

		    return selectedIndices;
		}


		Vector<float, DIM_X>& Vec_X() {return mVec_X;}
		const Vector<float, DIM_X>& Vec_X() const {return mVec_X;}
		const Vector<float, DIM_X> get_Vec_X() const {return mVec_X;}

		Matrix<float, DIM_X, DIM_X>& Mat_P() {return mMat_P;}
		const Matrix<float, DIM_X, DIM_X>& Mat_P() const {return mMat_P;}

		Matrix<float, DIM_X, DIM_X>& Mat_F() {return mMat_F;}
		const Matrix<float, DIM_X, DIM_X>& Mat_F() const {return mMat_F;}

		Matrix<float, DIM_X, DIM_X>& Mat_Q() {return mMat_Q;}
		const Matrix<float, DIM_X, DIM_X>& Mat_Q() const {return mMat_Q;}

		Matrix<float, DIM_Z, DIM_Z>& Mat_R() {return mMat_R;}
		const Matrix<float, DIM_Z, DIM_Z>& Mat_R() const {return mMat_R;}

		Matrix<float, DIM_Z, DIM_X>& Mat_H() {return mMat_H;}
		const Matrix<float, DIM_Z, DIM_X>& Mat_H() const {return mMat_H;}

		float Vec_X_arr[DIM_X][1];
		float Mat_P_arr[DIM_X][DIM_X];
		float Mat_F_arr[DIM_X][DIM_X];
		float Mat_Q_arr[DIM_X][DIM_X];
		float Mat_R_arr[DIM_Z][DIM_Z];
		float Mat_H_arr[DIM_Z][DIM_X];

	private:
		//State vector
		Vector<float, DIM_X> mVec_X;


		//State covariance matrix
		Matrix<float, DIM_X, DIM_X> mMat_P;


		//State transition matrix - A in neural decoding
		Matrix<float, DIM_X, DIM_X> mMat_F;


		//Process noise covariance matrix - W in neural decoding
		Matrix<float, DIM_X, DIM_X> mMat_Q;


		//Measurement noise covariance matrix - Q in neural decoding
		Matrix<float, DIM_Z, DIM_Z> mMat_R;


		//Measurement model matrix - H in neural decoding
		Matrix<float, DIM_Z, DIM_X> mMat_H;


		//Steady-state Kalman gain
		Matrix<float, DIM_X, DIM_Z> mMat_SSK;

		//Taylor approximation Kalman gain
		Matrix<float, DIM_X, DIM_Z> mMat_taylor_K;

		//Taylor approximation
		Matrix<float, DIM_Z, DIM_Z> mMat_C;

		float mInv_array1[DIM_Z][DIM_Z], mInv_array2[DIM_Z][DIM_Z];
		bool mPingpong = true;


		//For fitting
		MatrixXf mMat_Xtrain;
		MatrixXf mMat_Ytrain;

	};

}

#endif /* KALMANFILTER_H_ */

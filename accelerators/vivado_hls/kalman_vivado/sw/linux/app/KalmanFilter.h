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

#define N 10
#define TINY 1.0e-20
#define SIZE 1000


/* #include "A_array.h" */
/* #include "H_array.h" */
/* #include "W_array.h" */
/* #include "Q_array.h" */
/* #include "initial_state_array.h" */
/* #include "measurements_array.h" */
/* #include "prediction_array.h" */
/* #include "real_array.h" */


/* #include "A_array_soma.h" */
/* #include "H_array_soma.h" */
/* #include "W_array_soma.h" */
/* #include "Q_array_soma.h" */
/* #include "initial_state_array_soma.h" */
/* #include "measurements_array_soma.h" */
/* #include "prediction_array_soma.h" */
/* #include "real_array_soma.h" */

#include "A_array_hc.h"
#include "H_array_hc.h"
#include "W_array_hc.h"
#include "Q_array_hc.h"
#include "initial_state_array_hc.h"
#include "measurements_array_hc.h"
#include "prediction_array_hc.h"
#include "real_array_hc.h"


/* #include "inverse.h" */


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

		/* void correction_inv(const Vector<float, DIM_Z> &Vec_Z/\*measurement vector*\/ */
		/* 		//const Matrix<float, DIM_Z, DIM_X> &Mat_H/\*measurement model*\/, */
		/* 		//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/\*measurement noise covariance*\/ */
		/* 		) */
		/* { */
		/* 	//Innovation / residual */
		/* 	Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mVec_X); */

		/* 	//Innovation covariance */
		/* 	Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * mMat_P * mMat_H.transpose() + mMat_R; */

		/* 	//Invert the S matrix manually */
		/* 	Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_row = Mat_S; */
		/* 	//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_inv; */
		/* 	float S_array[DIM_Z][DIM_Z]; */
		/* 	float inv_array[DIM_Z][DIM_Z]; */
		/* 	Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row; */
		/* 	inverse_clean<DIM_Z>(S_array, inv_array); */

		/* 	Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&inv_array[0][0], DIM_Z, DIM_Z); */

		/* 	//Kalman gain */
		/* 	Matrix<float, DIM_X, DIM_Z> Mat_K = mMat_P * mMat_H.transpose() * Mat_S_inv; */

		/* 	//Identity matrix */
		/* 	Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity(); */

		/* 	//Update state vector */
		/* 	mVec_X = mVec_X  + Mat_K * Vec_Y; */

		/* 	//Update state covariance */
		/* 	mMat_P = (Mat_I -  Mat_K * mMat_H) * mMat_P; */

		/* } */

		/* void correction_inv_iter(const Vector<float, DIM_Z> &Vec_Z/\*measurement vector*\/, int iter */
		/* 		//const Matrix<float, DIM_Z, DIM_X> &Mat_H/\*measurement model*\/, */
		/* 		//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/\*measurement noise covariance*\/ */
		/* 		) */
		/* { */
		/* 	//Innovation / residual */
		/* 	Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mVec_X); */

		/* 	//Innovation covariance */
		/* 	Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * mMat_P * mMat_H.transpose() + mMat_R; */

		/* 	//Invert the S matrix manually */
		/* 	Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_row = Mat_S; */
		/* 	//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_inv; */
		/* 	float S_array[DIM_Z][DIM_Z]; */
		/* 	float inv_array[DIM_Z][DIM_Z], final_inv_array[DIM_Z][DIM_Z]; */
		/* 	Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row; */
		/* 	inverse_clean<DIM_Z>(S_array, inv_array); */

		/* 	many_iterative_inverse<DIM_Z>(S_array, inv_array, final_inv_array, iter); */

		/* 	Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&final_inv_array[0][0], DIM_Z, DIM_Z); */

		/* 	//Kalman gain */
		/* 	Matrix<float, DIM_X, DIM_Z> Mat_K = mMat_P * mMat_H.transpose() * Mat_S_inv; */

		/* 	//Identity matrix */
		/* 	Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity(); */

		/* 	//Update state vector */
		/* 	mVec_X = mVec_X  + Mat_K * Vec_Y; */

		/* 	//Update state covariance */
		/* 	mMat_P = (Mat_I -  Mat_K * mMat_H) * mMat_P; */

		/* } */

		/* void predict_and_correct_inv(const Vector<float, DIM_Z> &Vec_Z/\*measurement vector*\/ */
		/* 		//const Matrix<float, DIM_Z, DIM_X> &Mat_H/\*measurement model*\/, */
		/* 		//const Matrix<float, DIM_Z, DIM_Z> &Mat_R/\*measurement noise covariance*\/ */
		/* 		) */
		/* { */
		/* 	//Merge prediction into the computations of correction */

		/* 	//Innovation / residual */
		/* 	Vector<float, DIM_Z> Vec_Y = Vec_Z - (mMat_H * mMat_F * mVec_X); */

		/* 	//Innovation covariance */
		/* 	Matrix<float, DIM_Z, DIM_Z> Mat_S = mMat_H * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() + mMat_R; */

		/* 	//Invert the S matrix manually */
		/* 	Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_row = Mat_S; */
		/* 	//Matrix<float, DIM_Z, DIM_Z, RowMajor> Mat_S_inv; */
		/* 	float S_array[DIM_Z][DIM_Z]; */
		/* 	float inv_array[DIM_Z][DIM_Z]; */
		/* 	Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>>(&S_array[0][0], DIM_Z, DIM_Z) = Mat_S_row; */
		/* 	inverse_clean<DIM_Z>(S_array, inv_array); */

		/* 	Map<Matrix<float, DIM_Z, DIM_Z, RowMajor>> Mat_S_inv(&inv_array[0][0], DIM_Z, DIM_Z); */

		/* 	//Kalman gain */
		/* 	Matrix<float, DIM_X, DIM_Z> Mat_K = (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q) * mMat_H.transpose() * Mat_S_inv; */

		/* 	//Identity matrix */
		/* 	Matrix<float, DIM_X, DIM_X> Mat_I = Matrix<float, DIM_X, DIM_X>::Identity(); */

		/* 	//Update state vector */
		/* 	mVec_X = mMat_F * mVec_X  + Mat_K * Vec_Y; */

		/* 	//Update state covariance */
		/* 	mMat_P = (Mat_I -  Mat_K * mMat_H) * (mMat_F * mMat_P * mMat_F.transpose() + mMat_Q); */

		/* } */

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


		//For fitting
		MatrixXf mMat_Xtrain;
		MatrixXf mMat_Ytrain;

	};

}

#endif /* KALMANFILTER_H_ */

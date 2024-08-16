/*
 * clean_kalman.cpp
 *
OB *  Created on: Sep 17, 2023
 *      Author: geichler
 */

#include "c_run.h"
#include "KalmanFilter.h"

extern "C" {

void run_kalman()
{

	// int n = 0;
	// omp_set_num_threads(n);
	// Eigen::setNbThreads(n);
	// //Eigen::initParallel();
	// int nthreads = Eigen::nbThreads();

	std::cout << "\n\nTEST: Use Kalman filter on neural decoding data from the motor cortex\n\n" << std::endl;

	const int states = 6;
	const int neurons = 164;
	const int time_stamps = 10;//3793;
/*	const int states = 6;
	const int neurons = 52;
	const int time_stamps = 9199;*/
//	const int states = 6;
//	const int neurons = 46;
//	const int time_stamps = 4176;

	kf::KalmanFilterRegression<states/*state dimension*/, neurons/*neurons dimension*/> bci_kalman_filter;

	std::cout << "Initializing...\n";


	//Set initial position, velocity, acceleration
	Vector<float, states> vec_X = Map<Vector<float, states> >(initial);
	bci_kalman_filter.Vec_X() = vec_X;
//	for(int i = 0; i < states; i++)
//		bci_kalman_filter.Vec_X() = initial[i];

	std::cout << "X initialized\n";

	//Set initial state covariance matrix - zero matrix
	bci_kalman_filter.Mat_P() = Matrix<float, states, states>::Zero();

	std::cout << "P initialized\n";

	//Set state transition matrix
	Matrix<float, states, states> mat_F = Map<Matrix<float, states, states, RowMajor> >(A);
	bci_kalman_filter.Mat_F() = mat_F;
//	for(int i = 0; i < states*states; i++)
//		bci_kalman_filter.Mat_F() << A[i];

	std::cout << "F initialized\n";

	//Set noise covariance matrix
	Matrix<float, states, states> mat_Q = Map<Matrix<float, states, states, RowMajor> >(W);
	bci_kalman_filter.Mat_Q() = mat_Q;
//	for(int i = 0; i < states*states; i++)
//		bci_kalman_filter.Mat_Q() << W[i];

	std::cout << "Q initialized\n";

	//Set measurement noise covariance matrix
	Matrix<float, neurons, neurons> mat_R = Map<Matrix<float, neurons, neurons, RowMajor> >(Q);
	bci_kalman_filter.Mat_R() = mat_R;
//	for(int i = 0; i < neurons*neurons; i++)
//		bci_kalman_filter.Mat_R() << Q[i];

	std::cout << "R initialized\n";

	//Set measurement model matrix
	Matrix<float, neurons, states> mat_H = Map<Matrix<float, neurons, states, RowMajor> >(H);
	bci_kalman_filter.Mat_H() = mat_H;
//	for(int i = 0; i < neurons*states; i++)
//		bci_kalman_filter.Mat_H() << H[i];

	std::cout << "H initialized\n";


	float total_time_stamps = time_stamps;

	struct timespec startn, endn;
	double total_time = 0.0;


	float total_diff_inv = 0.0;
	//float total_time_stamps = 5/*time_stamps*/;
	float sum_r2_inv[6] = {0};
	float mean_inv[6] = {0};
	float sum_2_r2_inv[6] = {0};
	float R2[6] = {0};

	for(int j = 1; j < total_time_stamps/*time_stamps*/; j++)
	{

		Vector<float, neurons> vec_Z;
		for(int i = neurons*j; i < neurons*(j+1); i++)
			vec_Z(i-neurons*j) = measurements[i];

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startn);
		//bci_kalman_filter_inv.prediction();

		//bci_kalman_filter_inv.correction(vec_Z);
		//bci_kalman_filter_inv.correction_inv(vec_Z);
		//bci_kalman_filter_inv.correction_inv_iter(vec_Z, 1);


		bci_kalman_filter.predict_and_correct(vec_Z);
		//bci_kalman_filter_inv.predict_and_correct_opt(vec_Z);
		//bci_kalman_filter_inv.predict_and_correct_inv(vec_Z);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endn);

		float time_taken = (endn.tv_sec - startn.tv_sec) * 1e9;
		time_taken = (time_taken + (endn.tv_nsec - startn.tv_nsec)) * 1e-9;
		total_time += time_taken;

		//std::cout << "vector Z = \n" << vec_Z <<"\n";

//		std::cout << "After correction\n";
//		std::cout << "x = \n" << bci_kalman_filter.Vec_X() << "\n";
//		std::cout << "P = \n" << bci_kalman_filter.Mat_P() << "\n\n";

		Vector<float, states> vec_python; //first vector
		for(int i = states*j; i < states*(j+1); i++)
			vec_python(i-states*j) = prediction[i];

		Vector<float, states> vec_real_python; //first vector
		for(int i = states*j; i < states*(j+1); i++)
			vec_real_python(i-states*j) = real_out[i];

		Vector<float, states> diff_vec, diff_real;
		diff_vec = bci_kalman_filter.get_Vec_X() - vec_python;
		diff_real = bci_kalman_filter.get_Vec_X() - vec_real_python;
//		std::cout << "diff is = \n" << diff_vec << std::endl;
//		std::cout << "max is = \n" << diff_vec.maxCoeff() << std::endl;
//		std::cout << "min is = \n" << diff_vec.minCoeff() << std::endl;
		float max_diff = std::max(diff_vec.maxCoeff(), diff_vec.minCoeff()*(-1));
//		std::cout << "maximum diff = \n" << max_diff << std::endl;

//		std::cout << "After diff\n";
//		std::cout << "x = \n" << bci_kalman_filter_inv.Vec_X() << "\n";
//		std::cout << "python = \n" << vec_python << "\n";

		total_diff_inv += max_diff;

		for(int i = 0; i < states; i++){
			sum_r2_inv[i] += std::pow(diff_real(i), 2.0);
			mean_inv[i] += vec_real_python(i);
		}

		// std::cout << "finished iteration = \n" << j << "\n";

	}

	total_diff_inv = total_diff_inv / (total_time_stamps-1);
	std::cout << "Average diff is = \n" << total_diff_inv << std::endl;

	//Calculate R^2
	for(int i = 0; i < states; i++){
		mean_inv[i] = mean_inv[i] / (total_time_stamps-1);
	}

	for(int j = 1; j < total_time_stamps; j++){
		Vector<float, states> vec_real_python; //first vector
		for(int i = states*j; i < states*(j+1); i++)
			vec_real_python(i-states*j) = real_out[i];

		for(int i = 0; i < states; i++){
			sum_2_r2_inv[i] += std::pow(vec_real_python(i) - mean_inv[i], 2.0);
		}
	}

	//float R2[6] = {0};

	for(int i = 0; i < states; i++){
		R2[i] = 1-sum_r2_inv[i]/sum_2_r2_inv[i];
		printf("R2[%d] = %f\n", i, R2[i]);
	}

	printf("Total time for kalman filter: %f sec\n", total_time);

}

}  /* extern "C" */



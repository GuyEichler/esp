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
	const int time_stamps = 11;//3793;
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
        float max_diff_out = 0.0;
        float norm_max_diff_out = 0.0;
	//float total_time_stamps = 5/*time_stamps*/;
	float sum_r2_inv[6] = {0};
	float mean_inv[6] = {0};
	float sum_2_r2_inv[6] = {0};
	float R2[6] = {0};

        float sum_abs_vec = 0.0;
        float sum_sqr_vec = 0.0;

	for(int j = 1; j < total_time_stamps/*time_stamps*/; j++)
	{

		std::cout << "iteration start " << j <<"\n";

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

                Vector<float, states> diff_vec, diff_real, norm_diff_vec;
                diff_vec = bci_kalman_filter.get_Vec_X() - vec_python;
                for(int y = 0; y < states; y++)
                {
                    norm_diff_vec(y) =
                        std::abs(diff_vec(y)/(bci_kalman_filter.get_Vec_X())(y)) > std::abs(diff_vec(y)/vec_python(y)) ?
                        std::abs(diff_vec(y)/(bci_kalman_filter.get_Vec_X())(y)) : std::abs(diff_vec(y)/vec_python(y));
                }
                diff_real = bci_kalman_filter.get_Vec_X() - vec_real_python;

                float max_diff = std::max(diff_vec.maxCoeff(), diff_vec.minCoeff()*(-1));
                float norm_max_diff = std::max(norm_diff_vec.maxCoeff(), norm_diff_vec.minCoeff()*(-1));

                float abs_diff = (diff_vec.cwiseAbs()).sum();
                float sqr_diff = std::pow(abs_diff,2);
                sum_abs_vec += abs_diff;
                sum_sqr_vec += sqr_diff;

                total_diff_inv += max_diff;
                max_diff_out = std::max(max_diff_out, max_diff);
                norm_max_diff_out = std::max(norm_max_diff_out, norm_max_diff);

                for(int i = 0; i < states; i++){
                    sum_r2_inv[i] += std::pow(diff_real(i), 2.0);
                    mean_inv[i] += vec_real_python(i);
                }

		// std::cout << "finished iteration = \n" << j << "\n";
		std::cout << "iteration done " << j <<"\n";

	}

        sum_abs_vec = sum_abs_vec/((total_time_stamps-1)*states);
        sum_sqr_vec = sum_sqr_vec/((total_time_stamps-1)*states);

        total_diff_inv = total_diff_inv / (total_time_stamps-1);
        std::cout << "Average diff is = \n" << total_diff_inv << std::endl;
        std::cout << "Maximum diff is = \n" << max_diff_out << std::endl;
        std::cout << "Maximum norm diff is = \n" << norm_max_diff_out << std::endl;
        std::cout << "MAE is = \n" << sum_abs_vec << std::endl;
        std::cout << "MSE is = \n" << sum_sqr_vec << std::endl;

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

extern "C" {

void run_kalman_newton()
{

	// int n = 0;
	// omp_set_num_threads(n);
	// Eigen::setNbThreads(n);
	// //Eigen::initParallel();
	// int nthreads = Eigen::nbThreads();

	std::cout << "\n\nTEST: Use Kalman filter on neural decoding data from the motor cortex\n\n" << std::endl;

	const int states = 6;
	const int neurons = 164;
	const int time_stamps = 11;//3793;
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
        float max_diff_out = 0.0;
        float norm_max_diff_out = 0.0;
	//float total_time_stamps = 5/*time_stamps*/;
	float sum_r2_inv[6] = {0};
	float mean_inv[6] = {0};
	float sum_2_r2_inv[6] = {0};
	float R2[6] = {0};

        float sum_abs_vec = 0.0;
        float sum_sqr_vec = 0.0;

	for(int j = 1; j < total_time_stamps/*time_stamps*/; j++)
	{

		std::cout << "iteration start " << j <<"\n";

		Vector<float, neurons> vec_Z;
		for(int i = neurons*j; i < neurons*(j+1); i++)
			vec_Z(i-neurons*j) = measurements[i];

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startn);
		//bci_kalman_filter_inv.prediction();

		//bci_kalman_filter_inv.correction(vec_Z);
		//bci_kalman_filter_inv.correction_inv(vec_Z);
		//bci_kalman_filter_inv.correction_inv_iter(vec_Z, 1);


		//bci_kalman_filter.predict_and_correct(vec_Z);
		//bci_kalman_filter_inv.predict_and_correct_opt(vec_Z);
		//bci_kalman_filter_inv.predict_and_correct_inv(vec_Z);

                //Combinations with Newton - Gauss
                if(j == 1)
                    bci_kalman_filter.predict_and_correct_inv(vec_Z);
                else
                    bci_kalman_filter.predict_and_correct_iter(vec_Z, 1);


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

                Vector<float, states> diff_vec, diff_real, norm_diff_vec;
                diff_vec = bci_kalman_filter.get_Vec_X() - vec_python;
                for(int y = 0; y < states; y++)
                {
                    norm_diff_vec(y) =
                        std::abs(diff_vec(y)/(bci_kalman_filter.get_Vec_X())(y)) > std::abs(diff_vec(y)/vec_python(y)) ?
                        std::abs(diff_vec(y)/(bci_kalman_filter.get_Vec_X())(y)) : std::abs(diff_vec(y)/vec_python(y));
                }
                diff_real = bci_kalman_filter.get_Vec_X() - vec_real_python;

                float max_diff = std::max(diff_vec.maxCoeff(), diff_vec.minCoeff()*(-1));
                float norm_max_diff = std::max(norm_diff_vec.maxCoeff(), norm_diff_vec.minCoeff()*(-1));

                float abs_diff = (diff_vec.cwiseAbs()).sum();
                float sqr_diff = std::pow(abs_diff,2);
                sum_abs_vec += abs_diff;
                sum_sqr_vec += sqr_diff;

                total_diff_inv += max_diff;
                max_diff_out = std::max(max_diff_out, max_diff);
                norm_max_diff_out = std::max(norm_max_diff_out, norm_max_diff);

                for(int i = 0; i < states; i++){
                    sum_r2_inv[i] += std::pow(diff_real(i), 2.0);
                    mean_inv[i] += vec_real_python(i);
                }

		// std::cout << "finished iteration = \n" << j << "\n";
		std::cout << "iteration done " << j <<"\n";

	}

        sum_abs_vec = sum_abs_vec/((total_time_stamps-1)*states);
        sum_sqr_vec = sum_sqr_vec/((total_time_stamps-1)*states);

        total_diff_inv = total_diff_inv / (total_time_stamps-1);
        std::cout << "Average diff is = \n" << total_diff_inv << std::endl;
        std::cout << "Maximum diff is = \n" << max_diff_out << std::endl;
        std::cout << "Maximum norm diff is = \n" << norm_max_diff_out << std::endl;
        std::cout << "MAE is = \n" << sum_abs_vec << std::endl;
        std::cout << "MSE is = \n" << sum_sqr_vec << std::endl;

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

extern "C" {

void run_kalman_gauss()
{

	// int n = 0;
	// omp_set_num_threads(n);
	// Eigen::setNbThreads(n);
	// //Eigen::initParallel();
	// int nthreads = Eigen::nbThreads();

	std::cout << "\n\nTEST: Use Kalman filter on neural decoding data from the motor cortex\n\n" << std::endl;

	const int states = 6;
	const int neurons = 164;
	const int time_stamps = 11;//3793;
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
        float max_diff_out = 0.0;
        float norm_max_diff_out = 0.0;
	//float total_time_stamps = 5/*time_stamps*/;
	float sum_r2_inv[6] = {0};
	float mean_inv[6] = {0};
	float sum_2_r2_inv[6] = {0};
	float R2[6] = {0};

        float sum_abs_vec = 0.0;
        float sum_sqr_vec = 0.0;

	for(int j = 1; j < total_time_stamps/*time_stamps*/; j++)
	{

		std::cout << "iteration start " << j <<"\n";

		Vector<float, neurons> vec_Z;
		for(int i = neurons*j; i < neurons*(j+1); i++)
			vec_Z(i-neurons*j) = measurements[i];

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startn);
		//bci_kalman_filter_inv.prediction();

		//bci_kalman_filter_inv.correction(vec_Z);
		//bci_kalman_filter_inv.correction_inv(vec_Z);
		//bci_kalman_filter_inv.correction_inv_iter(vec_Z, 1);


		//bci_kalman_filter.predict_and_correct(vec_Z);
		//bci_kalman_filter_inv.predict_and_correct_opt(vec_Z);
		bci_kalman_filter.predict_and_correct_inv(vec_Z);
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

                Vector<float, states> diff_vec, diff_real, norm_diff_vec;
                diff_vec = bci_kalman_filter.get_Vec_X() - vec_python;
                for(int y = 0; y < states; y++)
                {
                    norm_diff_vec(y) =
                        std::abs(diff_vec(y)/(bci_kalman_filter.get_Vec_X())(y)) > std::abs(diff_vec(y)/vec_python(y)) ?
                        std::abs(diff_vec(y)/(bci_kalman_filter.get_Vec_X())(y)) : std::abs(diff_vec(y)/vec_python(y));
                }
                diff_real = bci_kalman_filter.get_Vec_X() - vec_real_python;

                float max_diff = std::max(diff_vec.maxCoeff(), diff_vec.minCoeff()*(-1));
                float norm_max_diff = std::max(norm_diff_vec.maxCoeff(), norm_diff_vec.minCoeff()*(-1));

                float abs_diff = (diff_vec.cwiseAbs()).sum();
                float sqr_diff = std::pow(abs_diff,2);
                sum_abs_vec += abs_diff;
                sum_sqr_vec += sqr_diff;

                total_diff_inv += max_diff;
                max_diff_out = std::max(max_diff_out, max_diff);
                norm_max_diff_out = std::max(norm_max_diff_out, norm_max_diff);

                for(int i = 0; i < states; i++){
                    sum_r2_inv[i] += std::pow(diff_real(i), 2.0);
                    mean_inv[i] += vec_real_python(i);
                }

		// std::cout << "finished iteration = \n" << j << "\n";
		std::cout << "iteration done " << j <<"\n";

	}

        sum_abs_vec = sum_abs_vec/((total_time_stamps-1)*states);
        sum_sqr_vec = sum_sqr_vec/((total_time_stamps-1)*states);

        total_diff_inv = total_diff_inv / (total_time_stamps-1);
        std::cout << "Average diff is = \n" << total_diff_inv << std::endl;
        std::cout << "Maximum diff is = \n" << max_diff_out << std::endl;
        std::cout << "Maximum norm diff is = \n" << norm_max_diff_out << std::endl;
        std::cout << "MAE is = \n" << sum_abs_vec << std::endl;
        std::cout << "MSE is = \n" << sum_sqr_vec << std::endl;

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


extern "C" {

void run_kalman_manual()
{

	// int n = 0;
	// omp_set_num_threads(n);
	// Eigen::setNbThreads(n);
	// //Eigen::initParallel();
	// int nthreads = Eigen::nbThreads();

	std::cout << "\n\nTEST: Use Kalman filter on neural decoding data from the motor cortex\n\n" << std::endl;

	const int states = 6;
	const int neurons = 164;
	const int time_stamps = 11;//3793;
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
        for(int i = 0; i < states; i++)
            bci_kalman_filter.Vec_X_arr[i][0] = vec_X(i);

	std::cout << "X initialized\n";

	//Set initial state covariance matrix - zero matrix
	bci_kalman_filter.Mat_P() = Matrix<float, states, states>::Zero();
        for(int i = 0; i < states; i++)
            for(int j = 0; j < states; j++)
                bci_kalman_filter.Mat_P_arr[i][j] = 0.0;

	std::cout << "P initialized\n";

	//Set state transition matrix
	Matrix<float, states, states> mat_F = Map<Matrix<float, states, states, RowMajor> >(A);
	bci_kalman_filter.Mat_F() = mat_F;
//	for(int i = 0; i < states*states; i++)
//		bci_kalman_filter.Mat_F() << A[i];
        for(int i = 0; i < states; i++)
            for(int j = 0; j < states; j++)
                bci_kalman_filter.Mat_F_arr[i][j] = mat_F(i,j);

	std::cout << "F initialized\n";

	//Set noise covariance matrix
	Matrix<float, states, states> mat_Q = Map<Matrix<float, states, states, RowMajor> >(W);
	bci_kalman_filter.Mat_Q() = mat_Q;
//	for(int i = 0; i < states*states; i++)
//		bci_kalman_filter.Mat_Q() << W[i];
        for(int i = 0; i < states; i++)
            for(int j = 0; j < states; j++)
                bci_kalman_filter.Mat_Q_arr[i][j] = mat_Q(i,j);

	std::cout << "Q initialized\n";

	//Set measurement noise covariance matrix
	Matrix<float, neurons, neurons> mat_R = Map<Matrix<float, neurons, neurons, RowMajor> >(Q);
	bci_kalman_filter.Mat_R() = mat_R;
//	for(int i = 0; i < neurons*neurons; i++)
//		bci_kalman_filter.Mat_R() << Q[i];
        for(int i = 0; i < neurons; i++)
            for(int j = 0; j < neurons; j++)
                bci_kalman_filter.Mat_R_arr[i][j] = mat_R(i,j);

	std::cout << "R initialized\n";

	//Set measurement model matrix
	Matrix<float, neurons, states> mat_H = Map<Matrix<float, neurons, states, RowMajor> >(H);
	bci_kalman_filter.Mat_H() = mat_H;
//	for(int i = 0; i < neurons*states; i++)
//		bci_kalman_filter.Mat_H() << H[i];
        for(int i = 0; i < neurons; i++)
            for(int j = 0; j < states; j++)
                bci_kalman_filter.Mat_H_arr[i][j] = mat_H(i,j);

	std::cout << "H initialized\n";


	float total_time_stamps = time_stamps;

	struct timespec startn, endn;
	double total_time = 0.0;


        float total_diff_inv = 0.0;
        float max_diff_out = 0.0;
        float norm_max_diff_out = 0.0;
	//float total_time_stamps = 5/*time_stamps*/;
	float sum_r2_inv[6] = {0};
	float mean_inv[6] = {0};
	float sum_2_r2_inv[6] = {0};
	float R2[6] = {0};

        float sum_abs_vec = 0.0;
        float sum_sqr_vec = 0.0;

	for(int j = 1; j < total_time_stamps/*time_stamps*/; j++)
	{

		std::cout << "iteration start " << j <<"\n";

		// Vector<float, neurons> vec_Z;
		// for(int i = neurons*j; i < neurons*(j+1); i++)
		// 	vec_Z(i-neurons*j) = measurements[i];

                float vec_Z_arr[neurons][1];
                for(int i = neurons*j; i < neurons*(j+1); i++)
                    vec_Z_arr[i-neurons*j][0] = measurements[i];

		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startn);
		//bci_kalman_filter_inv.prediction();

		//bci_kalman_filter_inv.correction(vec_Z);
		//bci_kalman_filter_inv.correction_inv(vec_Z);
		//bci_kalman_filter_inv.correction_inv_iter(vec_Z, 1);


		//bci_kalman_filter.predict_and_correct(vec_Z);
		//bci_kalman_filter_inv.predict_and_correct_opt(vec_Z);
		//bci_kalman_filter.predict_and_correct_inv(vec_Z);
                bci_kalman_filter.predict_and_correct_inv_manual(vec_Z_arr);
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

                Vector<float, states> diff_vec, diff_real, norm_diff_vec;
                Map<Vector<float, states>> X_diff(&bci_kalman_filter.Vec_X_arr[0][0]);
                diff_vec = X_diff - vec_python;
                // diff_vec = bci_kalman_filter.get_Vec_X() - vec_python;
                for(int y = 0; y < states; y++)
                {
                    norm_diff_vec(y) =
                        std::abs(diff_vec(y)/X_diff(y)) > std::abs(diff_vec(y)/vec_python(y)) ?
                        std::abs(diff_vec(y)/X_diff(y)) : std::abs(diff_vec(y)/vec_python(y));
                }
                diff_real = X_diff - vec_real_python;

                float max_diff = std::max(diff_vec.maxCoeff(), diff_vec.minCoeff()*(-1));
                float norm_max_diff = std::max(norm_diff_vec.maxCoeff(), norm_diff_vec.minCoeff()*(-1));

                float abs_diff = (diff_vec.cwiseAbs()).sum();
                float sqr_diff = std::pow(abs_diff,2);
                sum_abs_vec += abs_diff;
                sum_sqr_vec += sqr_diff;

                total_diff_inv += max_diff;
                max_diff_out = std::max(max_diff_out, max_diff);
                norm_max_diff_out = std::max(norm_max_diff_out, norm_max_diff);

                for(int i = 0; i < states; i++){
                    sum_r2_inv[i] += std::pow(diff_real(i), 2.0);
                    mean_inv[i] += vec_real_python(i);
                }

		// std::cout << "finished iteration = \n" << j << "\n";
		std::cout << "iteration done " << j <<"\n";

	}

        sum_abs_vec = sum_abs_vec/((total_time_stamps-1)*states);
        sum_sqr_vec = sum_sqr_vec/((total_time_stamps-1)*states);

        total_diff_inv = total_diff_inv / (total_time_stamps-1);
        std::cout << "Average diff is = \n" << total_diff_inv << std::endl;
        std::cout << "Maximum diff is = \n" << max_diff_out << std::endl;
        std::cout << "Maximum norm diff is = \n" << norm_max_diff_out << std::endl;
        std::cout << "MAE is = \n" << sum_abs_vec << std::endl;
        std::cout << "MSE is = \n" << sum_sqr_vec << std::endl;

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

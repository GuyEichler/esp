#include "c_run.h"
#include "eigen/Eigen/Dense"
// #include "eigen/Eigen/SVD"
#include "gemm_stratus.h"

// #include "cfg.h"
// #include <libesp.h>

#include <time.h>
#include <bits/stdc++.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

/* <<--params-def-->> */
#define DO_RELU    0
#define TRANSPOSE  1
#define NINPUTS    2
#define D3         8
#define D2         8
#define D1         8
#define ST_OFFSET  (NINPUTS * (D1 * D2 + D2 * D3))
#define LD_OFFSET1 0
#define LD_OFFSET2 (NINPUTS * (D1 * D2))

extern "C" void esp_dummy_gemm(void *x);

using namespace Eigen;
using namespace std;

// -- from utility.h and utility.cc of ekf
namespace utility
{
enum SensorType { LASER, RADAR };

struct SensorReading {
    long long       timestamp;
    SensorType      sensor_type;
    Eigen::VectorXd measurement;
};

const Eigen::VectorXd CalculateRmse(const std::vector<Eigen::VectorXd> &estimations,
                                    const std::vector<Eigen::VectorXd> &ground_truth);
const Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd &x_state);
const Eigen::VectorXd PolarToCartesian(const Eigen::VectorXd &polar_vector);
const Eigen::VectorXd CartesianToPolar(const Eigen::VectorXd &x_state);

void CheckArguments(int argc, char *argv[]);
void CheckFiles(std::ifstream &in_file, std::string &in_name, std::ofstream &out_file, std::string &out_name);
}; // namespace utility

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace utility
{
void CheckArguments(int argc, char *argv[])
{
    std::string usage_instructions = "Usage instructions: ";
    usage_instructions += argv[0];
    usage_instructions += " path/to/input.txt output.txt";

    bool has_valid_args = false;

    // make sure the user has provided input and output files
    if (argc == 1) {
        std::cerr << usage_instructions << std::endl;
    } else if (argc == 2) {
        std::cerr << "Please include an output file.\n" << usage_instructions << std::endl;
    } else if (argc == 3) {
        has_valid_args = true;
    } else if (argc > 3) {
        std::cerr << "Too many arguments.\n" << usage_instructions << std::endl;
    }

    if (!has_valid_args) {
        exit(EXIT_FAILURE);
    }
}

void CheckFiles(ifstream &in_file, string &in_name, ofstream &out_file, string &out_name)
{
    if (!in_file.is_open()) {
        std::cerr << "Cannot open input file: " << in_name << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!out_file.is_open()) {
        std::cerr << "Cannot open output file: " << out_name << std::endl;
        exit(EXIT_FAILURE);
    }
}

const VectorXd CalculateRmse(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth)
{

    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;
    if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
        cerr << "The estimation vector size should equal ground truth\
                        vector size. Also, the estimation vector size should\
                        not be zero"
             << endl;
        return rmse;
    }

    for (int i = 0; i < estimations.size(); ++i) {
        VectorXd residual = estimations[i] - ground_truth[i];
        residual          = residual.array() * residual.array();

        rmse += residual;
    }

    rmse /= estimations.size();
    rmse = rmse.array().sqrt();

    return rmse;
}

const MatrixXd CalculateJacobian(const VectorXd &x_state)
{

    MatrixXd Hj = MatrixXd::Zero(3, 4);

    // recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    // pre-compute a set of terms to avoid repeated calculation
    float c1 = px * px + py * py;
    float c2 = sqrt(c1);
    float c3 = (c1 * c2);

    // check division by zero
    if (fabs(c1) < 1e-4) {
        return Hj;
    }

    // compute the Jacobian matrix
    Hj << (px / c2), (py / c2), 0, 0, -(py / c1), (px / c1), 0, 0, py * (vx * py - vy * px) / c3,
        px * (px * vy - py * vx) / c3, px / c2, py / c2;

    return Hj;
}

const VectorXd PolarToCartesian(const VectorXd &polar_vector)
{
    VectorXd cart_vector(4);

    float rho     = polar_vector(0);
    float phi     = polar_vector(1);
    float rho_dot = polar_vector(2);

    float x  = rho * cos(phi);
    float y  = rho * sin(phi);
    float vx = rho_dot * cos(phi);
    float vy = rho_dot * sin(phi);

    cart_vector << x, y, vx, vy;
    return cart_vector;
}

const Eigen::VectorXd CartesianToPolar(const Eigen::VectorXd &x_state)
{
    VectorXd polar_vector(3);

    float x  = x_state(0);
    float y  = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    if (fabs(x + y) < 1e-4) {
        x = 1e-4;
        y = 1e-4;
    }

    float rho     = sqrt(x * x + y * y);
    float phi     = atan2(y, x);
    float rho_dot = (x * vx + y * vy) / rho;

    polar_vector << rho, phi, rho_dot;
    return polar_vector;
}

}; // namespace utility
// -- end of utility.h and utility.cc of ekf

void *    cfg;
unsigned *do_relu_i;
unsigned *transpose_i;
unsigned *ninputs_i;
unsigned *d3_i;
unsigned *d2_i;
unsigned *d1_i;
unsigned *st_offset_i;
unsigned *ld_offset1_i;
unsigned *ld_offset2_i;
unsigned *src_offset_i;
unsigned *dst_offset_i;
int *     acc_buf_i;

#define FX_IL 16

static inline float fixed32_to_float(int value, int n_int_bits)
{
    unsigned shift_int = 0x3f800000 - 0x800000 * (32 - n_int_bits);
    float *  shift     = (float *)(&shift_int);

    return (*shift) * (float)value;
}

static inline int float_to_fixed32(float value, int n_int_bits)
{
    unsigned shift_int = 0x3f800000 + 0x800000 * (32 - n_int_bits);
    float *  shift     = (float *)&shift_int;

    return (int)(value * (*shift));
}

template <typename MatrixType> MatrixType operator*(const MatrixType &A, const MatrixType &B)
{
    // cout << "CUSTOM" << endl;

    // esp_run(cfg_000, 1);

    int rowsA = A.rows();
    int colsA = A.cols();
    int rowsB = B.rows();
    int colsB = B.cols();

    // Check if the matrices can be multiplied
    if (colsA != rowsB) {
        std::cerr << "Error: The number of columns in A must match the number of rows in B." << std::endl;
        return MatrixType();
    }

    *do_relu_i   = 0;
    *transpose_i = 0;
    *ninputs_i   = 1;

    *d1_i = rowsA;
    *d2_i = colsA;
    *d3_i = colsB;

    *ld_offset1_i = 0;
    *ld_offset2_i = 1 * rowsA * colsA;
    *st_offset_i  = rowsA * colsA + rowsB * colsB;

    *src_offset_i = 0;
    *dst_offset_i = 0;

    // Allocate memory for the arrays
    int *A_data_int = &(acc_buf_i[0]);               // new int[rowsA * colsA];
    int *B_data_int = &(acc_buf_i[(*ld_offset2_i)]); // new int[rowsB * colsB];
    int *C_data_int = &(acc_buf_i[(*st_offset_i)]);  // new int[rowsA * colsB];

    // Copy the data from the matrices into the arrays
    for (int i = 0; i < rowsA; ++i)
        for (int j = 0; j < colsA; ++j) {
            A_data_int[i * colsA + j] = float_to_fixed32(A(i, j), FX_IL);
        }
    for (int i = 0; i < rowsB; ++i)
        for (int j = 0; j < colsB; ++j) {
            B_data_int[i * colsB + j] = float_to_fixed32(B(i, j), FX_IL);
        }

    // acc_buf_i = acc_buf;

    // Run gemm accelerator
    esp_dummy_gemm(cfg);

    // Copy the result from the array into a matrix
    MatrixType C_acc(rowsA, colsB);
    for (int i = 0; i < rowsA; ++i)
        for (int j = 0; j < colsB; ++j) {
            C_acc(i, j) = fixed32_to_float(C_data_int[i * colsB + j], FX_IL);
        }

    return C_acc;
}

template <typename MatrixType> MatrixType CustomProductC(const MatrixType &A, const MatrixType &B)
{
    // cout << "CUSTOM" << endl;

    // esp_run(cfg_000, 1);

    int rowsA = A.rows();
    int colsA = A.cols();
    int rowsB = B.rows();
    int colsB = B.cols();

    // Check if the matrices can be multiplied
    if (colsA != rowsB) {
        std::cerr << "Error: The number of columns in A must match the number of rows in B." << std::endl;
        return MatrixType();
    }

    // Allocate memory for the arrays
    double *A_data = new double[rowsA * colsA];
    double *B_data = new double[rowsB * colsB];
    double *C_data = new double[rowsA * colsB];

    // Copy the data from the matrices into the arrays
    for (int i = 0; i < rowsA; ++i)
        for (int j = 0; j < colsA; ++j) {
            A_data[i * colsA + j] = A(i, j);
        }
    for (int i = 0; i < rowsB; ++i)
        for (int j = 0; j < colsB; ++j) {
            B_data[i * colsB + j] = B(i, j);
        }

    // Perform the matrix multiplication using a custom algorithm
    for (int i = 0; i < rowsA; ++i)
        for (int j = 0; j < colsB; ++j) {
            double sum = 0.0;
            for (int k = 0; k < colsA; ++k)
                sum += A_data[i * colsA + k] * B_data[k * colsB + j];
            C_data[i * colsB + j] = sum;
        }

    // Copy the result from the array into a matrix
    MatrixType C(rowsA, colsB);
    for (int i = 0; i < rowsA; ++i)
        for (int j = 0; j < colsB; ++j) {
            C(i, j) = C_data[i * colsB + j];
        }

    // Deallocate memory for the arrays
    delete[] A_data;
    delete[] B_data;
    delete[] C_data;

    return C;
}

// extern "C" {

void c_run_gemm(int m, int n, int p, void *x, unsigned *do_relu, unsigned *transpose, unsigned *ninputs, unsigned *d1, unsigned *d2,
                unsigned *d3, unsigned *st_offset, unsigned *ld_offset1, unsigned *ld_offset2, unsigned *src_offset,
                unsigned *dst_offset, int *acc_buf)
{

    // set global pointers to gemm accelerator configuration
    cfg = x;

    do_relu_i    = do_relu;
    transpose_i  = transpose;
    ninputs_i    = ninputs;
    d3_i         = d3;
    d2_i         = d2;
    d1_i         = d1;
    st_offset_i  = st_offset;
    ld_offset1_i = ld_offset1;
    ld_offset2_i = ld_offset2;
    src_offset_i = src_offset;
    dst_offset_i = dst_offset;
    acc_buf_i    = acc_buf;

    // MatrixXf A(2, 2);
    // A << 1, 2,
    //     3, 4;

    // MatrixXf B(2, 2);
    // B << 5, 6,
    //     7, 8;

    // MatrixXf(row, col)
    cout << "(m, n, p) = (" << m << ", " << n << ", " << p << ")" << endl;
    MatrixXf A = MatrixXf::Random(m, n);
    MatrixXf B = MatrixXf::Random(n, p);

    // cout << "A : " << A << endl;
    // cout << "B : " << B << endl;

    cout << setprecision(15) << endl;

    struct timespec startn, endn;
    struct timespec startn1, endn1;
    struct timespec startn2, endn2;

    // -- [Test 0]: run gemm accelerator
    // cout << "[Test 0]: Using the custom implementation and overloaded operator (gemm accelerator): " << endl;

    clock_gettime(CLOCK_MONOTONIC, &startn);
    MatrixXf C0 = A * B;
    clock_gettime(CLOCK_MONOTONIC, &endn);
    // cout << "The product of A and B is:\n" << C0 << endl;

    double time_taken = (endn.tv_sec - startn.tv_sec) * 1e9;
    time_taken        = (time_taken + (endn.tv_nsec - startn.tv_nsec)) * 1e-9;
    cout << fixed << "[Test 0] Time taken: " << time_taken << " sec" << endl;

    // -- [Test 1]: run on riscv
    // cout << "[Test 1]: Using the custom implementation in regular C: " << endl;

    clock_gettime(CLOCK_MONOTONIC, &startn1);
    MatrixXf C1 = CustomProductC(A, B);
    clock_gettime(CLOCK_MONOTONIC, &endn1);
    // cout << "The product of A and B is:\n" << C1 << endl;

    double time_taken1 = (endn1.tv_sec - startn1.tv_sec) * 1e9;
    time_taken1        = (time_taken1 + (endn1.tv_nsec - startn1.tv_nsec)) * 1e-9;
    cout << fixed << "[Test 1] Time taken: " << time_taken1 << " sec" << endl;

    // -- [Test 2]: run on riscv with Eigen implementation
    // cout << "[Test 2]: Using the Eigen implementation and original operator: " << endl;

    clock_gettime(CLOCK_MONOTONIC, &startn2);
    MatrixXf C2 = A.operator*(B);
    clock_gettime(CLOCK_MONOTONIC, &endn2);
    // cout << "The product of A and B is:\n" << C2 << endl;

    double time_taken2 = (endn2.tv_sec - startn2.tv_sec) * 1e9;
    time_taken2        = (time_taken2 + (endn2.tv_nsec - startn2.tv_nsec)) * 1e-9;
    cout << fixed << "[Test 2] Time taken: " << time_taken2 << " sec" << endl;

    // Eigen::MatrixXf C2 = Eigen::MatrixBase<Eigen::MatrixXf>::operator*(A,B);

    MatrixXf C01 = C0 - C1;

    // cout << endl;
    // cout << "Maximum difference between accelerator and C: " << C01.maxCoeff() << endl;
    // cout << "Minimum difference between accelerator and C: " << C01.minCoeff() << endl;

    // cout << "== c_run_gemm done ==" << endl;
    // cout << endl;
}

// -- from ekf.h and ekf.cc of ekf
class ExtendedKalmanFilter
{
  public:
    // state vector
    Eigen::VectorXd x_;

    // state covariance matrix
    Eigen::MatrixXd P_;

    // state transistion matrix
    Eigen::MatrixXd F_;

    // process covariance matrix
    Eigen::MatrixXd Q_;

    // measurement matrix
    Eigen::MatrixXd H_;

    // measurement covariance matrix
    Eigen::MatrixXd R_;

    // 4x4 Identity matrix we will need later.
    Eigen::MatrixXd I_;

    /**
     * Prediction Predicts the state and the state covariance
     * using the process model
     * @param delta_T Time between k and k+1 in s
     */
    void Predict();

    /**
     * Updates the state by using standard Kalman Filter equations
     * @param z The measurement at k+1
     */
    void Update(const Eigen::VectorXd &z);
    void UpdateEkf(const Eigen::VectorXd &z);

  private:
    // Keep things DRY
    void CallRestOfUpdate(const Eigen::VectorXd &z);
};

using utility::CartesianToPolar;

void ExtendedKalmanFilter::Predict()
{
    x_ = F_ * x_;
    // [kuanlin]: matrix multiplication here:
    P_ = F_ * P_ * F_.transpose() + Q_;
    // P_ = (F_.operator*(P_)).operator*(F_.transpose()) + Q_;
}

void ExtendedKalmanFilter::Update(const VectorXd &z)
{
    VectorXd y = z - H_ * x_; // error calculation
    CallRestOfUpdate(y);
}

void ExtendedKalmanFilter::UpdateEkf(const VectorXd &z)
{
    VectorXd hx = CartesianToPolar(x_);
    VectorXd y  = z - hx;
    CallRestOfUpdate(y);
}

void ExtendedKalmanFilter::CallRestOfUpdate(const VectorXd &y)
{
    MatrixXd Ht = H_.transpose();

    // [kuanlin]: matrix multiplication here:
    MatrixXd PHt = P_ * Ht;
    MatrixXd S   = H_ * PHt + R_;
    // MatrixXd PHt = P_.operator*(Ht);
    // MatrixXd S   = H_.operator*(PHt) + R_;

    MatrixXd K = PHt * S.inverse();

    // New state
    x_ = x_ + (K * y);

    // [kuanlin]: matrix multiplication here:
    P_ = (I_ - K * H_) * P_;
    // P_ = (I_ - K.operator*(H_)).operator*(P_);
}
// -- end of ekf.h and ekf.cc of ekf

// -- from fusion_ekf.h and fusion_ekf.cc of ekf
class FusionEkf
{

  public:
    FusionEkf();
    void            ProcessMeasurement(utility::SensorReading &reading);
    Eigen::VectorXd current_estimate();

  private:
    bool      is_initialized_;
    long long previous_timestamp_;

    // This Fusion class handles LIDAR and RADAR
    // So we set the sensor measurement maxtrices and covar for the sensors
    Eigen::MatrixXd R_laser_;
    Eigen::MatrixXd H_laser_;
    Eigen::MatrixXd R_radar_;

    // since this is a constant velocity model. we add the accel as noise.
    float noise_ax_;
    float noise_ay_;

    /**
     * Extended Kalman Filter update and prediction math lives in here.
     */
    ExtendedKalmanFilter ekf_;

    void UpdateFQ(float dt);
};

using utility::SensorReading;
using utility::SensorType;

FusionEkf::FusionEkf()
{
    is_initialized_     = false;
    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);

    // measurement covariance matrix - laser
    R_laser_ << 0.0225, 0, 0, 0.0225;

    // measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;

    H_laser_ << 1, 0, 0, 0, 0, 1, 0, 0;

    noise_ax_ = 9;
    noise_ay_ = 9;

    // setup efk F
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1;

    // setup efk P
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 1000;

    ekf_.Q_ = MatrixXd::Zero(4, 4);
    ekf_.I_ = MatrixXd::Identity(4, 4);
}

Eigen::VectorXd FusionEkf::current_estimate() { return ekf_.x_; }

void FusionEkf::UpdateFQ(float dt)
{
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    // Modify the F matrix so that the time is integrated
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    // Modify the Q matrix so that the time is integrated
    ekf_.Q_(0, 0) = (dt_4 / 4) * noise_ax_;
    ekf_.Q_(0, 2) = (dt_3 / 2) * noise_ax_;

    ekf_.Q_(1, 1) = (dt_4 / 4) * noise_ay_;
    ekf_.Q_(1, 3) = (dt_3 / 2) * noise_ay_;

    ekf_.Q_(2, 0) = (dt_3 / 2) * noise_ax_;
    ekf_.Q_(2, 2) = dt_2 * noise_ax_;

    ekf_.Q_(3, 1) = (dt_3 / 2) * noise_ay_;
    ekf_.Q_(3, 3) = dt_2 * noise_ay_;
}

void FusionEkf::ProcessMeasurement(SensorReading &reading)
{
    if (!is_initialized_) {
        // we will need to init the state vector x,y,vx,vy
        switch (reading.sensor_type) {
            case SensorType::LASER:
                ekf_.x_ = VectorXd(4);
                ekf_.x_ << reading.measurement[0], reading.measurement[1], 0, 0;
                break;
            case SensorType::RADAR:
                ekf_.x_ = utility::PolarToCartesian(reading.measurement);
                break;
        }

        // we can't let x and y be zero.
        if (fabs(ekf_.x_(0) + ekf_.x_(1)) < 1e-4) {
            ekf_.x_(0) = 1e-4;
            ekf_.x_(1) = 1e-4;
        }

        previous_timestamp_ = reading.timestamp;
        is_initialized_     = true;
        return;
    }

    float dt            = (reading.timestamp - previous_timestamp_) / 1000000.0; // dt - expressed in seconds
    previous_timestamp_ = reading.timestamp;

    // We update the F and Q matrices using dt
    UpdateFQ(dt);

    // Make Prediction only if we the dt is big enough.
    if (dt > 1e-3) {
        ekf_.Predict();
    }

    switch (reading.sensor_type) {
        case SensorType::RADAR:
            ekf_.H_ = utility::CalculateJacobian(ekf_.x_);
            ekf_.R_ = R_radar_;
            ekf_.UpdateEkf(reading.measurement);
            break;

        case SensorType::LASER:
            ekf_.H_ = H_laser_;
            ekf_.R_ = R_laser_;
            ekf_.Update(reading.measurement);
            break;
    }
} // end FusionEkf::ProcessMeasurement

// -- end of fusion_ekf.h and fusion_ekf.cc of ekf

using utility::CalculateRmse;
using utility::CheckArguments;
using utility::CheckFiles;
using utility::SensorReading;
using utility::SensorType;

void c_run_ekf(int argc, char **argv, void *x, unsigned *do_relu, unsigned *transpose, unsigned *ninputs, unsigned *d1,
               unsigned *d2, unsigned *d3, unsigned *st_offset, unsigned *ld_offset1, unsigned *ld_offset2,
               unsigned *src_offset, unsigned *dst_offset, int *acc_buf)
{
    struct timespec start_ekf, end_ekf;


    CheckArguments(argc, argv);

    string   in_file_name_ = argv[1];
    ifstream in_file_(in_file_name_.c_str(), ifstream::in);

    string   out_file_name_ = argv[2];
    ofstream out_file_(out_file_name_.c_str(), ofstream::out);

    CheckFiles(in_file_, in_file_name_, out_file_, out_file_name_);

    vector<SensorReading> sensor_readings;
    vector<VectorXd>      ground_truths;

    string line;

    while (getline(in_file_, line)) {
        istringstream iss(line);

        string        sensor_type;
        SensorReading sensor_reading;
        long long     timestamp;

        // reads first element from the current line
        iss >> sensor_type;
        if (sensor_type.compare("L") == 0) {
            // LASER MEASUREMENT
            float x, y;
            iss >> x;
            iss >> y;
            iss >> timestamp;

            sensor_reading.sensor_type = SensorType::LASER;
            sensor_reading.measurement = VectorXd(2);
            sensor_reading.measurement << x, y;
            sensor_reading.timestamp = timestamp;
        } else if (sensor_type.compare("R") == 0) {
            // RADAR MEASUREMENT
            float ro, phi, ro_dot;
            iss >> ro;
            iss >> phi;
            iss >> ro_dot;
            iss >> timestamp;

            sensor_reading.sensor_type = SensorType::RADAR;
            sensor_reading.measurement = VectorXd(3);
            sensor_reading.measurement << ro, phi, ro_dot;
            sensor_reading.timestamp = timestamp;
        }

        sensor_readings.push_back(sensor_reading);
        // read ground truth data to compare later
        float x_gt, y_gt, vx_gt, vy_gt;
        iss >> x_gt;
        iss >> y_gt;
        iss >> vx_gt;
        iss >> vy_gt;

        VectorXd ground_truth(4);
        ground_truth << x_gt, y_gt, vx_gt, vy_gt;
        ground_truths.push_back(ground_truth);
    }

    FusionEkf        fusion_ekf_;
    vector<VectorXd> estimations;


    clock_gettime(CLOCK_MONOTONIC, &start_ekf);


    size_t N = sensor_readings.size();
    for (size_t k = 0; k < N; ++k) {
        SensorReading reading = sensor_readings[k];
        fusion_ekf_.ProcessMeasurement(reading);

        // output the estimation
        VectorXd x_ = fusion_ekf_.current_estimate();

        out_file_ << x_(0) << "\t";
        out_file_ << x_(1) << "\t";
        out_file_ << x_(2) << "\t";
        out_file_ << x_(3) << "\t";

        switch (reading.sensor_type) {
            case SensorType::RADAR:
                out_file_ << reading.measurement(0) * cos(reading.measurement(1)) << "\t";
                out_file_ << reading.measurement(0) * sin(reading.measurement(1)) << "\t";
                break;

            case SensorType::LASER:
                out_file_ << reading.measurement(0) << "\t";
                out_file_ << reading.measurement(1) << "\t";
                break;
        }

        // output the ground truth packages
        VectorXd ground_truth = ground_truths[k];
        out_file_ << ground_truth(0) << "\t";
        out_file_ << ground_truth(1) << "\t";
        out_file_ << ground_truth(2) << "\t";
        out_file_ << ground_truth(3) << "\n";

        estimations.push_back(x_);
    }

    // compute the accuracy (RMSE)
    cout << "Accuracy - RMSE:" << endl << CalculateRmse(estimations, ground_truths) << endl;

    clock_gettime(CLOCK_MONOTONIC, &end_ekf);
 double time_taken = (end_ekf.tv_sec - start_ekf.tv_sec) * 1e9;
    time_taken        = (time_taken + (end_ekf.tv_nsec - start_ekf.tv_nsec)) * 1e-9;
    cout << fixed << "Time taken: " << time_taken << " sec\n" << endl;


    // close files
    if (out_file_.is_open()) {
        out_file_.close();
    }

    if (in_file_.is_open()) {
        in_file_.close();
    }
}

// } /* extern "C" */

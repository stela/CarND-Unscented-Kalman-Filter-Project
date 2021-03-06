#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Array2d;
using Eigen::Array3d;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() :
    // Process noise standard deviation longitudinal acceleration in m/s^2
    // best found so far: RMSE x, y, vx, vy
    // a=1.5, yawdd=0.5:  0.0692, 0.0829, 0.3333, 0.2345
    // Max RMSE allowed: .09, .10, .40, .30

    std_a_(1.5),

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_(0.5),

    //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
    // Laser measurement noise standard deviation position1 in m
    std_laspx_(0.15),

    // Laser measurement noise standard deviation position2 in m
    std_laspy_(0.15),

    // Radar measurement noise standard deviation radius in m
    std_radr_(0.3),

    // Radar measurement noise standard deviation angle in rad
    std_radphi_(0.03),

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_(0.3),
    //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

    // State dimension (vector length)
    n_x_(5),

    // Augmented state dimension
    n_aug_(n_x_ + 2),

    // Number of sigma points
    n_sig_(2 * n_aug_ + 1),

    lambda_(3 - n_x_)
  {

  // weights of sigma points
  weights_ = sigmaPointWeights(n_x_, n_aug_, n_sig_, lambda_);

  // initial state vector, should be overwritten by first measurement
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  Xsig_pred_ = MatrixXd::Zero(n_x_, n_sig_);

  R_radar_ = MatrixXd::Zero(3, 3);
  Array3d r_diag(std_radr_ * std_radr_, std_radphi_ * std_radphi_, std_radrd_ * std_radrd_);
  R_radar_.diagonal() = r_diag;

  R_lidar_ = MatrixXd::Zero(2, 2);
  Array2d l_diag(std_laspx_ * std_laspx_, std_laspy_ * std_laspy_);
  R_lidar_.diagonal() = l_diag;
}

VectorXd UKF::sigmaPointWeights(const int n_x, const int n_aug, const int n_sig, const double lambda) const {
  //create vector for weights
  VectorXd weights = VectorXd(n_sig);

  // set weights
  double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<n_sig; i++) {
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }
  return weights;
}


UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  //
  // Initialization similar to EKF project
  //

  if (!is_initialized_) {
    // first measurement
    x_ = VectorXd(n_x_);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = measurement_pack.raw_measurements_(0);     // distance
      double theta = measurement_pack.raw_measurements_(1);   // bearing/angle
      double rho_dot = measurement_pack.raw_measurements_(2); // radial velocity (delta-ro)
      long long ts = measurement_pack.timestamp_;
      x_(0) = rho * cos(theta);  // x
      x_(1) = rho * sin(theta);  // y
      x_(2) = rho_dot;            // velocity
      x_(3) = 0.0;  // yaw angle unknown
      x_(4) = 0.0;  // yaw angle rate/velocity unknown

      cout << "UKF: first measurement is RADAR" << endl << x_ << endl;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      double px = measurement_pack.raw_measurements_(0);
      double py = measurement_pack.raw_measurements_(1);
      long long ts = measurement_pack.timestamp_;
      // Data set 1 has a Laser line first, contains these ground_truth values:
      // x_groundtruth, y_groundtruth, vx_groundtruth, vy_groundtruth, yaw_groundtruth, yawrate_groundtruth
      // 6.000000e-01,  6.000000e-01,  5.199937e+00,   0,              0,               6.911322e-03
      // delta_t appears to be 100000 us (=0.1s) to the next (radar) measurement
      // the RMSE will be very high for the first few measurements, but let's ignore those ;-)
      x_(0) = px;
      x_(1) = py;
      x_(2) = 0.0; // radial velocity unknown
      x_(3) = 0.0; // yaw angle unknown
      x_(4) = 0.0; // yaw angle rate/velocity unknown

      cout << "UKF: first measurement is LASER" << endl << x_ << endl;
    } else {
      cout << "ERROR! Unexpected first measurement type: " << measurement_pack.sensor_type_ << endl;
    }

    time_us_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  } // end !is_initialized


  //
  // Control structure similar to EKF project - Predict then Update
  //

  // Prediction
  double dt = (measurement_pack.timestamp_ - time_us_) / 1000000.0;
  Prediction(dt);

  // Update
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(measurement_pack);
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
    UpdateLidar(measurement_pack);
  } else {
    cout << "UKF: WARNING! Unknown sensor " << measurement_pack.sensor_type_ << endl;
  }
  time_us_ = measurement_pack.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  //
  // AUGMENTED SIGMA POINT GENERATION
  //

  // Sigma point generator from "Lesson 7: 15. Generating Sigma Points Assignment 2",
  // then expanded upon in Lesson 7: 18. Augmentation Assignment 2"

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  //set first column of sigma point matrix
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  //
  // SIGMA POINT PREDICTION
  //

  // From "Lesson 7: 21. Sigma Point Prediction Assignment 2":

  //predict sigma points

  for (int i = 0; i< n_sig_; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }


  //
  // PREDICT MEAN AND COVARIANCE
  //

  // From "Lesson 7: 24. Predicted mean and covariance assignment 2"

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3));

    P_ += weights_(i) * x_diff * x_diff.transpose() ;
  }
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  // Basic KF, no extended/UKF required (just like previous EKF project)
  // See "Lesson 6: 13. Laser Measurements Part 4"
  // but instead of H being 2x4, it's now 2x5 (append zeroes)
  // Print out state + covariance to debug

  // laser measurement is 2D, store in vector for easier processing
  const int n_z = 2;

  // transform the sigma points into measurement space
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, n_sig_);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_sig_; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  S += R_lidar_;

  // 2. Update state

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  for (int i = 0; i < n_sig_; i++) {  // 2n+1 sigma points

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // Incoming radar measurement
  VectorXd z = meas_package.raw_measurements_;
  // residual
  VectorXd z_diff = z - z_pred;
  NormalizeAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K*S*K.transpose();

  // NIS Lidar Update
  NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  //
  // PREDICT RADAR SIGMA POINTS
  //

  // From "Lesson 7: 27. Predict Radar Measurement Assignment 2"

  //set measurement dimension, radar can measure r, phi, and r_dot
  const int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  //transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++) {  //2n+1 sigma points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  for (int i=0; i < n_sig_; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  for (int i = 0; i < n_sig_; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
      0, std_radphi_*std_radphi_, 0,
      0, 0,std_radrd_*std_radrd_;
  S += R;


  //
  // UKF Update (radar)
  //

  // From "Lesson 7: 30. UKF Update Assignment 2"

  VectorXd z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3));

    Tc += + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;
  NormalizeAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K*S*K.transpose();
}

// Normalize angle in-place, make sure reference to actual location and not temporary storage
void UKF::NormalizeAngle(double &phi) {
  // alternative, but suspect slower and less numerically stable, for small angle corrections at least
  // phi = atan2(sin(phi), cos(phi));

  while (phi >  M_PI) phi -= 2.*M_PI;
  while (phi < -M_PI) phi += 2.*M_PI;
}

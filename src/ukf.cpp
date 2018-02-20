#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() :
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_(30),

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_(30),

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

    // weights of sigma points
    weights_(sigmaPointWeights())
  {
  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.setOnes();



  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

VectorXd UKF::sigmaPointWeights() {
  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
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
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  //
  // Initialization similar to EKF project
  //

  if (!is_initialized_) {
    // first measurement
    x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double ro = measurement_pack.raw_measurements_(0);
      double theta = measurement_pack.raw_measurements_(1);
      double ro_dot = measurement_pack.raw_measurements_(2);
      long long ts = measurement_pack.timestamp_;
      x_(0) = ro * cos(theta);
      x_(1) = ro * sin(theta);
      x_(2) = ro_dot * cos(theta);
      x_(3) = ro_dot * sin(theta);

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
      x_(2) = 0.0;
      x_(3) = 0.0;

      cout << "UKF: first measurement is LASER" << endl << x_ << endl;
    }

    // TODO assign to this.time_us_ ???
    time_us_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  } // end !is_initialized


  //
  // Control structure similar to EKF project
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

  // print the output
  //cout << "x_ = " << x_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //
  // AUGMENTED SIGMA POINT GENERATION
  //

  // Sigma point generator from "Lesson 7: 15. Generating Sigma Points Assignment 2",
  // then expanded upon in Lesson 7: 18. Augmentation Assignment 2"

  //define spreading parameter
  const double lambda = 3 - n_x_;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  //set first column of sigma point matrix
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda+n_aug_) * L.col(i);
  }

  //
  // SIGMA POINT PREDICTION
  //

  // From "Lesson 7: 21. Sigma Point Prediction Assignment 2":
  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //predict sigma points

  for (int i = 0; i< 2*n_aug_+1; i++)
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
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }


  //
  // PREDICT MEAN AND COVARIANCE
  //

  // From "Lesson 7: 24. Predicted mean and covariance assignment 2"

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);

  // set weights
  double weight_0 = lambda/(lambda+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda);
    weights(i) = weight;
  }

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    // TODO change to P += to avoid IDE warning
    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }

  //
  // ASSIGN OUTPUTS
  //

  // TODO add missing field updates
  Xsig_pred_ = Xsig_pred;
  x_ = x;
  P_ = P;
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Makes use of lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // Basic KF, no extended/UKF required (just like previous EKF project)
  // See "Lesson 6: 13. Laser Measurements Part 4"
  // but instead of H being 2x4, it's now 2x5 (append zeroes)
  // Print out state + covariance to debug
  // Do normalize all the angles (already done?)

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

    // TODO normalize angle z_diff(1) here?

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  S = S + R_lidar_;

  // 2. Update state
  // TODO move z decl down to first use
  // Incoming radar measurement
  VectorXd z = meas_package.raw_measurements_;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // TODO normalize angle z_diff(1) here?

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // TODO normalize angle __x_diff(3)__ here?

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;
  // TODO normalize angle z_diff(1) here?

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // NIS Lidar Update
  NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //
  // PREDICT RADAR SIGMA POINTS
  //

  // From "Lesson 7: 27. Predict Radar Measurement Assignment 2"

  //set measurement dimension, radar can measure r, phi, and r_dot
  const int n_z = 3;

  //radar measurement noise standard deviation radius in m
  const double std_radr = 0.3;

  //radar measurement noise standard deviation angle in rad
  const double std_radphi = 0.0175;

  //radar measurement noise standard deviation radius change in m/s
  const double std_radrd = 0.1;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

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
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr*std_radr, 0, 0,
      0, std_radphi*std_radphi, 0,
      0, 0,std_radrd*std_radrd;
  S = S + R;


  //
  // UKF Update (radar)
  //

  // From "Lesson 7: 30. UKF Update Assignment 2"

  //define spreading parameter
  double lambda = 3 - n_aug_;

  VectorXd z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  // TODO replace 2*n_aug_+1 with n_sig, everywhere
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    // TODO extract angle normalization into reusable function
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}

// Normalize angle in-place, make sure reference to actual location and not temporary storage
void UKF::NormalizeAngle(double &phi) {
  // alternative, suspect slower and less numerically stable, for small angle corrections at least
  // phi = atan2(sin(phi), cos(phi));

  while (phi >  M_PI) phi -= 2.*M_PI;
  while (phi < -M_PI) phi += 2.*M_PI;
}

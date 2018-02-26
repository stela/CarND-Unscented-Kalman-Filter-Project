#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_ = false;

  ///* if this is false, laser measurements will be ignored (except for init)
  const bool use_laser_ = true;

  ///* if this is false, radar measurements will be ignored (except for init)
  const bool use_radar_ = true;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  // uses the Constant Turn Rat and Velocity Magnitude Model (CTRV) model for state
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  // Equivalent to previous_timestamp_ in the EKF project
  long long time_us_ = 0;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  const double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  const double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  const double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  const double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  const double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  const double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  const double std_radrd_ ;

  ///* Weights of sigma points
  const VectorXd weights_;

  ///* State dimension
  const int n_x_;

  ///* Augmented state dimension
  const int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

private:
  ///* Number of sigma points
  const int n_sig_;

  ///* Noise covariance matrix of radar measurements
  MatrixXd R_radar_;

  ///* Noise covariance matrix of lidar measurements
  MatrixXd R_lidar_;

  ///* NIS for radar
  double NIS_radar_;

  ///* NIS for lidar
  double NIS_lidar_;





public:

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

private:
  /**
   * Normalizes the angle given to between +- M_PI
   * @param phi reference to angle to normalize
   */
  void NormalizeAngle(double &phi);

  /**
   * Creates weights of sigma points
   */
   VectorXd sigmaPointWeights();
};

#endif /* UKF_H */

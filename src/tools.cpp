#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {
  max_rmse = VectorXd(4);
  max_rmse << 0,0,0,0;
}

Tools::~Tools() = default

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // From "Lesson 6: 23. Ec=valuating KF Performance 2"
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.empty()) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse / estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // update the max components
  bool updated = false;
  for (int i=0; i < rmse.size(); i++) {
    if (rmse(i) > max_rmse(i)) {
      max_rmse(i) = rmse(i);
      updated = true;
    }
  }
  if (updated) {
    cout << "New RMSE high! " << endl << max_rmse << endl;
    cout << "Current RMSE: " << endl << rmse << endl;
  }

  // return the result
  return rmse;
 }

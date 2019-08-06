#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::string;
using std::vector;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft_ = F_.transpose();
  P_ = F_ * P_ * Ft_ + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * update the state by using Kalman Filter equations
   */
  
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * update the state by using Extended Kalman Filter equations
   */
  
  VectorXd z_pred(3);
  float x2y2 = (x_[0] * x_[0]) + (x_[1] * x_[1]);  // px² + py²
  // transform from carthesian to polar coodinate system / h(x')
  if (x2y2 != 0){
      z_pred << sqrt(x2y2), // rho
  				atan2(x_[1], x_[0]), // phi
  				((x_[0]*x_[2]) + (x_[1]*x_[3])) / sqrt(x2y2); //rho dot
  }
  else {
      std::cout << "carth -> polar - Error - Div/0" << std::endl;
      z_pred << sqrt(x2y2), // rho
  				atan2(x_[1], x_[0]), // phi
  				((x_[0]*x_[2]) + (x_[1]*x_[3])) / 0.1; //rho dot
  }
  
  VectorXd y = z - z_pred;
  
  // normalize phi
  float pi = 3.14159265359;
  while (y[1] > pi){
    std::cout << "normalize phi" << std::endl;
    y[1] -= pi*2;
  }
  while (y[1] < (pi*-1)){
    std::cout << "normalize phi" << std::endl;
    y[1] += pi*2;
  }
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
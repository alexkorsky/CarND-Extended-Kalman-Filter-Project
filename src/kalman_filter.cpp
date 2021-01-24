#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  // direct from the Project Quizes:

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // direct from the Project Quizes:

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
   * TODO: update the state by using Extended Kalman Filter equations
   */


	/*

	In the radar update step, the Jacobian matrix H_j
	is used to calculate S, K and P. To calculate y, we use the equations that map the predicted location x'
	from Cartesian coordinates to polar coordinates.

	The predicted measurement vector x'is a vector containing values in the form  p_x, p_y, v_x, v_y
	The radar sensor will output values in polar coordinates: rho, phi, rho_dot

	In order to calculate y for the radar sensor, we need to convert x'to polar coordinates.
	In other words, the function h(x) maps values from Cartesian coordinates to polar coordinates.
	So the equation for radar becomes y = z_radar - h(x')

	*/

	  float rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
	  float phi = atan2(x_(1), x_(0));
	  float rho_dot;

	  if (fabs(rho) < 0.0001) {//avoid devision by zero
	    rho_dot = 0;
	  } else {
	    rho_dot = (x_(0)*x_(2) + x_(1)*x_(3))/rho;
	  }

	  VectorXd h_x(3);
	  h_x << rho, phi, rho_dot;

	  //then use the same KF equations
	  VectorXd y = z - h_x;
	  MatrixXd Ht = H_.transpose();
	  MatrixXd S = H_ * P_ * Ht + R_;
	  MatrixXd Si = S.inverse();
	  MatrixXd PHt = P_ * Ht;
	  MatrixXd K = PHt * Si;

	  //normalizing angles till in [-pi : pi] range
	  while (y(1) > M_PI) y(1) -= 2*M_PI;
	  while (y(1) < -M_PI) y(1) += 2*M_PI;


	  //new estimate
	  x_ = x_ + (K * y);
	  long x_size = x_.size();
	  MatrixXd I = MatrixXd::Identity(x_size, x_size);
	  P_ = (I - K * H_) * P_;
}

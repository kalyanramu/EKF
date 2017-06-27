#include "kalman_filter.h"
#include <math.h>
#include <iostream>
#define _USE_MATH_DEFINES 

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
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
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];
  
  cout << "X vector: " << endl << x_ << endl;
  float rho = sqrt(px*px+py*py);
  // avoid division by zero
  float theta;
  if(fabs(px) < 0.001){
    cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << endl;
  }
  else {
  theta = atan2(py, px);
  }
  //float theta = atan2(py,px); //*Need to convert the angle to -pi,pi
  cout << "theta:" << theta << endl;
  float rho_dot;
  if (rho < 0.001)
  {
    cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << endl;
  }
  else
  {
    rho_dot = ((px*vx) + (py*vy))/rho;
  }
  
  VectorXd zpred(3);
  zpred << rho,theta,rho_dot;
  VectorXd y = z- zpred; //can exceed pi (e.g., pi -(-pi))


//Normalize error angle  
  while (y(1)> M_PI) {
    y(1) -= 2 * M_PI;   
  }
  
  while (y(1)<-M_PI) {
    y(1) += 2 * M_PI;
  }
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  //new estimate/state
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(x_.size(),x_.size());//https://eigen.tuxfamily.org/dox/group__TutorialAdvancedInitialization.html
  P_ = (I - K * H_) * P_;
  
}








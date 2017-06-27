#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

#define EPS 0.0001 // A very small number
#define EPS2 0.0000001

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  //VectorXd diff(4) ;
  	for(int i=0; i < estimations.size(); ++i){
        // ... your code here
  		VectorXd residual = estimations[i] - ground_truth[i];
  		residual = residual.array()*residual.array();
  		rmse += residual;
		
	}

	rmse = rmse/estimations.size(); //Divide by number of measurements

	//calculate the square root
	rmse = rmse.array().sqrt();
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
    // Code from lectures quizes
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  MatrixXd Hj(3,4);
  // Deal with the special case problems
  if (fabs(px) < EPS and fabs(py) < EPS){
	  px = EPS;
	  py = EPS;
  }
  // Pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;
  // Check division by zero
  if(fabs(c1) < EPS){
  cout << "Calculate Jacobian  -- Error division by zero" << endl;
  //getchar();
	c1 = EPS;
  }
  float c2 = sqrt(c1);
  float c3 = (c1*c2);
  // Compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
       -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  return Hj;
}

#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#define EPS 0.0001 // A very small number

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
//Note: ekf_ is the instance of kalmanfilter inside FusionEKF instance (refer to FusionEKF.h)

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  long long previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
 cout << "Fusion EKF Initialized" << endl;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 0, 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float r0 = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      float rdot = measurement_pack.raw_measurements_[2];
      ekf_.x_ << r0*cos(theta), r0*sin(theta),rdot*cos(theta),rdot*sin(theta) ; //* Need to figure our velocity
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];
      ekf_.x_ << px,py,0,0;
    }
     // Deal with the special case initialisation problems
    if (fabs(ekf_.x_(0)) < EPS and fabs(ekf_.x_(1)) < EPS){
		ekf_.x_(0) = EPS;
		ekf_.x_(1) = EPS;
		}

    //Initialize State Transition Matrix
    ekf_.F_ = MatrixXd(4,4);
    ekf_.F_ << 1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1;

    //Initialize Co-variance Matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
			   0, 1, 0, 0,
			   0, 0, 500, 0,
			   0, 0, 0, 500;

		//Initialize Process Co-varance Matrix
		ekf_.Q_ = MatrixXd::Zero(4,4);


    previous_timestamp_ = measurement_pack.timestamp_; //update timestamp
    // done initializing, no need to predict or update
    is_initialized_ = true;
    cout << "EKF Initialiation Done";
    
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  //Update F, Q matrices before calling Predict [To do] --- Use Lesson 5, Laser Measurements Part 3 equations
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_; //update timestamp
  //cout << "dt: "<<dt << endl;
  float dt2 = dt*dt;
  float dt3 = dt2*dt;
  float dt4 = dt3*dt;
  //update F
  ekf_.F_ << 1,0,dt,0,
  						0,1,0,dt,
  						0,0,1,0,
  						0,0,0,1;

  //update Q
  float noise_ax = 9.0; float noise_ay = 9.0;
  ekf_.Q_ << (dt4*noise_ax/4.0),0,(dt3*noise_ax/2.0),0,
  						0, (dt4*noise_ay/4.0),0,(dt3*noise_ay/2.0),
  						(dt3*noise_ax/2.0),0,(dt2*noise_ax),0,
  						0,(dt3*noise_ay/2.0),0,(dt2*noise_ay);


  //cout << "Before EKF Predicted";
  ekf_.Predict();
  //cout << "After EKF Predicted";

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_); //use calaculate jacobian here [To do]
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    
    
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

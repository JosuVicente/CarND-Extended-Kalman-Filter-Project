#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
* Constructor.
*/
FusionEKF::FusionEKF() {
	is_initialized_ = false;

	previous_timestamp_ = 0;

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

	H_laser_ << 1, 0, 0, 0,
		0, 1, 0, 0;

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
		* Initialize the state ekf_.x_ with the first measurement.
		**/

		// first measurement
		cout << "EKF: " << endl;
		ekf_.x_ = VectorXd(4);
		ekf_.x_ << 1, 1, 1, 1;

		//first timestamp
		previous_timestamp_ = measurement_pack.timestamp_;

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			float rho = measurement_pack.raw_measurements_[0];
			float phi = measurement_pack.raw_measurements_[1];
			float rhodot = measurement_pack.raw_measurements_[2];

			float px = rho * cos(phi);
			float py = rho * sin(phi);
			float vx = rhodot * cos(phi);
			float vy = rhodot * sin(phi);

			if (px == 0 && py == 0) {
				return;
			}

			ekf_.x_ << px, py, vx, vy;
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			*/
			float px = measurement_pack.raw_measurements_[0];
			float py = measurement_pack.raw_measurements_[1];

			if (px == 0 && py == 0) {
				return;
			}

			ekf_.x_ << px, py, 0, 0;
		}

		/**
		* Create the state covariance matrix P.
		*/
		ekf_.P_ = MatrixXd(4, 4);
		ekf_.P_ << 1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1000, 0,
			0, 0, 0, 1000;

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

	/*****************************************************************************
	*  Prediction
	****************************************************************************/

	/**
	Compute the time elapsed between the current and previous measurements */
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;
	/* Update the state transition matrix F according to the new elapsed time .*/
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, dt, 0,
		0, 1, 0, dt,
		0, 0, 1, 0,
		0, 0, 0, 1;

	/* Update the process noise covariance matrix.*/
	float noise_ax = 9;
	float noise_ay = 9;

	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << pow(dt, 4)*noise_ax / 4, 0, pow(dt, 3)*noise_ax / 2, 0,
		0, pow(dt, 4)*noise_ay / 4, 0, pow(dt, 3)*noise_ay / 2,
		pow(dt, 3)*noise_ax / 2, 0, pow(dt, 2)*noise_ax, 0,
		0, pow(dt, 3)*noise_ay / 2, 0, pow(dt, 2)*noise_ay;


	ekf_.Predict();

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
		Hj_ = tools.CalculateJacobian(ekf_.x_);
		ekf_.H_ = Hj_;
		ekf_.R_ = R_radar_;
		ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	}
	else {
		// Laser updates	  
		ekf_.H_ = H_laser_;
		ekf_.R_ = R_laser_;
		ekf_.Update(measurement_pack.raw_measurements_);
	}

	// print the output
	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}

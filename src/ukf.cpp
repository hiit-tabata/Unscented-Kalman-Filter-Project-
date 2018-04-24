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
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;

  n_aug_ = 7;

  n_radar_ = 3;

  n_laser_ = 2;

  lambda_ = 3 - n_aug_;

  x_aug_ = VectorXd(n_aug_);

  deltax_ = VectorXd(n_aug_);

  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  P_aug_ = MatrixXd(n_aug_, n_aug_);

  L_ = MatrixXd(n_aug_, n_aug_);

  // init weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++)
    weights_(i) = 0.5 / (n_aug_ + lambda_);

  // Variables used for radar update
  z_pred_radar_ = VectorXd(n_radar_);
  deltaz_radar_ = VectorXd(n_radar_);
  Zsig_radar_ = MatrixXd(n_radar_, 2 * n_aug_ + 1);
  // Radar measurement noise covariance matrix is constant/persistent
  R_radar_ = MatrixXd(n_radar_, n_radar_);
  R_radar_.fill(0.);
  R_radar_(0, 0) = std_radr_ * std_radr_;
  R_radar_(1, 1) = std_radphi_ * std_radphi_;
  R_radar_(2, 2) = std_radrd_ * std_radrd_;
  S_radar_ = MatrixXd(n_radar_, n_radar_);
  Tc_radar_ = MatrixXd(n_x_, n_radar_);
  K_radar_ = MatrixXd(n_x_, n_radar_);

  // Variables used for laser update
  z_pred_laser_ = VectorXd(n_laser_);
  deltaz_laser_ = VectorXd(n_laser_);
  Zsig_laser_ = MatrixXd(n_laser_, 2 * n_aug_ + 1);
  // Laser measurement noise covariance matrix is constant/persistent
  R_laser_ = MatrixXd(n_laser_, n_laser_);
  R_laser_.fill(0.);
  R_laser_(0, 0) = std_laspx_ * std_laspx_;
  R_laser_(1, 1) = std_laspy_ * std_laspy_;
  S_laser_ = MatrixXd(n_laser_, n_laser_);
  Tc_laser_ = MatrixXd(n_x_, n_laser_);
  K_laser_ = MatrixXd(n_x_, n_laser_);

  // State covariance matrix P will be initialized using some data
  // from the first measurement, in ProcessMeasurement() below.
}

UKF::~UKF()
{
  // close the nis reader
  NISvals_radar_.close();
  NISvals_laser_.close();
}


void UKF::initialize(MeasurementPackage meas_package)
{
  // init the state of ekf_.x with the first measurement.
  // create the covariance matrix
  cout << "init. unscented kalman Filter" << endl;

  // Handle the radar init.
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    double rho = meas_package.raw_measurements_[0];
    double phi = meas_package.raw_measurements_[1];

    x_ << rho * cos(phi), rho * sin(phi), 0., 0., 0.;
  }
  else if( meas_package.sensor_type_ == MeasurementPackage::LASER ) 
  {
    // handle laser case
    x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0., 0., 0.;
  }
  previous_timestamp_ = meas_package.timestamp_;

  // init the state covariance matrix
  P_.fill(0.);
  P_(0, 0) = 1.;
  P_(1, 1) = 1.;
  P_(2, 2) = 1.;
  P_(3, 3) = 1.;
  P_(4, 4) = 1.;
  
  is_initialized_ = true; 
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  
  if (!is_initialized_)
  {
    initialize(meas_package);
    return;
  }

  double dt = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  if(dt > 0.0001){
    Prediction(dt);
  }

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    UpdateRadar(meas_package);
  }else{
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{

  // create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  // create augmented covariance matrix
  P_aug_.fill(0.);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;

  // create square root matrix
  L_ = P_aug_.llt().matrixL();

  // create augmented sigma points 
  Xsig_aug_.col(0) = x_aug_;
  for(int i = 0; i < n_aug_; i++){
    Xsig_aug_.col(1 +i) = x_aug_ + sqrt(lambda_ + n_aug_) * L_.col(i);
    Xsig_aug_.col(1 +i + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L_.col(i);
  }

  // augmented sigma points created
  // loop all augmented sigma points to process nodel with noise 
  // do the prediction steps
  for (int pt = 0; pt <= 2*n_aug_; pt++){      
      //extract values for better readability
      double p_x      = Xsig_aug_(0,pt);
      double p_y      = Xsig_aug_(1,pt);
      double v        = Xsig_aug_(2,pt);
      double yaw      = Xsig_aug_(3,pt);
      double yawd     = Xsig_aug_(4,pt);
      double nu_a     = Xsig_aug_(5,pt);
      double nu_yawdd = Xsig_aug_(6,pt);
      
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
      Xsig_pred_(0,pt) = px_p;
      Xsig_pred_(1,pt) = py_p;
      Xsig_pred_(2,pt) = v_p;
      Xsig_pred_(3,pt) = yaw_p;
      Xsig_pred_(4,pt) = yawd_p;
  }


  // weight init. in contructor 

  // predict state mean 
  x_.fill(0.);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // predict state covariance matrix
  P_.fill(0.);
  for (int i = 0; i <= 2 * n_aug_; i ++){
    deltax_ = Xsig_pred_.col(i) - x_;
    deltax_ = Xsig_pred_.col(i) - x_;
    while( deltax_(3) > M_PI ){
      deltax_(3) -= 2.*M_PI;
    }
    while( deltax_(3) < -M_PI ){
       deltax_(3) += 2.*M_PI;
    }

    P_ = P_ + weights_(i)*deltax_*deltax_.transpose();  
  }
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_pack
 */
// Much of this code is taken from previously completed quizzes.
void UKF::UpdateRadar(MeasurementPackage meas_pack) 
{
  // Transform sigma points into measurement space
  // 2n+1 simga points
  for( int pt = 0; pt < 2*n_aug_ + 1; pt++ )
  {
    // extract values for better readibility
    double p_x  = Xsig_pred_(0,pt);
    double p_y  = Xsig_pred_(1,pt);
    double v    = Xsig_pred_(2,pt);
    double yaw  = Xsig_pred_(3,pt);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_radar_(0,pt) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_radar_(1,pt) = atan2(p_y,p_x);                                 //phi
    Zsig_radar_(2,pt) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }
  
  //mean predicted measurement
  z_pred_radar_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_radar_ = z_pred_radar_ + weights_(i) * Zsig_radar_.col(i);
  }
  
  while( z_pred_radar_(1) > M_PI ) z_pred_radar_(1)-=2.*M_PI;
  while( z_pred_radar_(1) <-M_PI ) z_pred_radar_(1)+=2.*M_PI;

  S_radar_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S_radar_ = S_radar_ + weights_(i) * z_diff * z_diff.transpose();
  }

  // got the S fir radar

  S_radar_ = S_radar_ + R_radar_;

  //calculate cross correlation matrix
  Tc_radar_.fill(0.);
  for( int pt = 0; pt < 2*n_aug_ + 1; pt++ )
  {
      deltax_ = Xsig_pred_.col(pt) - x_;
      deltaz_radar_ = Zsig_radar_.col(pt) - z_pred_radar_;
      while( deltax_(1)> M_PI ) deltax_(1)-=2.*M_PI;
      while( deltax_(1)<-M_PI ) deltax_(1)+=2.*M_PI;
      while( deltaz_radar_(1)> M_PI ) deltaz_radar_(1)-=2.*M_PI;
      while( deltaz_radar_(1)<-M_PI ) deltaz_radar_(1)+=2.*M_PI;
      Tc_radar_ = Tc_radar_ + weights_(pt)*deltax_*deltaz_radar_.transpose();
  }
  
  //calculate K_radar_alman gain K_radar_;
  K_radar_ = Tc_radar_*S_radar_.inverse();
  
  //update state mean and covariance matrix
  deltaz_radar_ = meas_pack.raw_measurements_ - z_pred_radar_;
  while( deltaz_radar_(1) > M_PI ) deltaz_radar_(1)-=2.*M_PI;
  while( deltaz_radar_(1) <-M_PI ) deltaz_radar_(1)+=2.*M_PI;
  
  x_ = x_ + K_radar_*deltaz_radar_;
  P_ = P_ - K_radar_*S_radar_*K_radar_.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_pack
 */
void UKF::UpdateLidar(MeasurementPackage meas_pack) 
{
  // Transform sigma points into measurement space
  for( int pt = 0; pt < 2*n_aug_ + 1; pt++ )
  {
    Zsig_laser_(0,pt) = Xsig_pred_(0,pt);
    Zsig_laser_(1,pt) = Xsig_pred_(1,pt);
  }
  
  //calculate mean predicted measurement
  z_pred_laser_.fill(0.);
  for( int pt = 0; pt < 2*n_aug_ + 1; pt++ ){
    z_pred_laser_ = z_pred_laser_ + weights_(pt)*Zsig_laser_.col(pt);
  }
  
  //calculate measurement covariance matrix S_laser_
  S_laser_.fill(0.);
  for( int pt = 0; pt < 2*n_aug_ + 1; pt++ )
  {
    deltaz_laser_ = Zsig_laser_.col(pt) - z_pred_laser_;
    S_laser_ = S_laser_ + weights_(pt)*deltaz_laser_*deltaz_laser_.transpose();
  }

  S_laser_ = S_laser_ + R_laser_;

  //calculate cross correlation matrix
  Tc_laser_.fill(0.);
  for( int pt = 0; pt < 2*n_aug_ + 1; pt++ )
  {
      deltax_ = Xsig_pred_.col(pt) - x_;
      deltaz_laser_ = Zsig_laser_.col(pt) - z_pred_laser_;
      Tc_laser_ = Tc_laser_ + weights_(pt)*deltax_*deltaz_laser_.transpose();
  }
  
  //calculate K_laser_alman gain K_laser_;
  K_laser_ = Tc_laser_*S_laser_.inverse();
  
  //update state mean and covariance matrix
  deltaz_laser_ = meas_pack.raw_measurements_ - z_pred_laser_;
  
  // update state
  x_ = x_ + K_laser_*deltaz_laser_;
  // update covariance matrix
  P_ = P_ - K_laser_*S_laser_*K_laser_.transpose();
}

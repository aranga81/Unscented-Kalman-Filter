#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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

  // flag to check if state estimation via UKF initialized
  is_initialized_ = false;

  // number of states in CTRV model 
  n_x_ = 5;

  // Augmented states - [n_x_ nu_a nu_yawdd] - 7 for CTRV
  n_aug_ = n_x_ + 2;

  // lambda scaling parameter for sigma point generation
  lambda_ = 3 - n_x_;

  // time stamp initialization - microseconds
  time_us_ = 0.0;

  // weights for scaling back sigma point to mean & covariance
  weights_ = VectorXd(2 * n_aug_ + 1);

  //NIS 
  NIS_laser_ = 0.0;

  NIS_ = 0.0;

  NIS_radar_ = 0.0;

  // Predicted state sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /*************************************************************
  *                   INITIALIZATION      
  *************************************************************/

  if(!is_initialized_)
  {

    time_us_ = meas_package.timestamp_;

    x_ << 1., 1., 1., 0., 0.;

    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    {
      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);

      float x = rho * cos(phi);
      float y = rho * sin(phi);
      float vx = 0.0;
      float vy = 0.0;
      float v = 0.0;
      float psi = 0.0;
      float psi_dot = 0.0;

      x_ << x, y, v, psi, psi_dot; 

    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    {

      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0.0, 0.0, 0.0;

    }

    is_initialized_ = true;

    cout << "UKF initial state: " << x_ <<endl; 
    return;
  }

  /*************************************************************
  *                   PREDICTION      
  *************************************************************/

  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0; 
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);


  /*************************************************************
  *                   UPDATE      
  *************************************************************/

  if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    UpdateLidar(meas_package);
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }

  cout << "UKF state: " << x_ <<endl;
  cout << "UKF covariance: " << P_ <<endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  /*************************************************************
  *                   Generate Sigma Points       
  *************************************************************/

  MatrixXd X_sigma = MatrixXd(n_x_, 2 * n_x_ + 1);
  lambda_ = 3 - n_x_;

  //Calculate Square root of Covariance P - Cholesky decomposition
  MatrixXd A = P_.llt().matrixL();

  X_sigma.col(0) = x_;

  for(int i=0; i<n_x_; i++)
  {
    X_sigma.col(i + 1) = x_ + sqrt(lambda_ + n_x_) * A.col(i);
    X_sigma.col(n_x_ + i + 1) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
  }
  
  /*************************************************************
  *       Augment Sigma Points with process noise      
  *************************************************************/
  VectorXd x_aug = VectorXd(n_aug_);

  MatrixXd X_sigma_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  lambda_ = 3 - n_aug_;

  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // Augmented Covariance Matrix P
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //Calculate Square root of Covariance P_aug - Cholesky decomposition
  MatrixXd M = P_aug.llt().matrixL();

  //Generate Augmented Sigma Points
  X_sigma_aug.col(0) = x_aug;

  for(int i=0; i<n_aug_; i++)
  {
    X_sigma_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * M.col(i);
    X_sigma_aug.col(n_aug_+i+1) = x_aug - sqrt(lambda_ + n_aug_) * M.col(i);
  }

  /*************************************************************
  *           Prediction Sigma Points     
  *************************************************************/

  for (int i = 0; i < 2*n_aug_+1; i++)
  {
    double px = X_sigma_aug(0, i);
    double py = X_sigma_aug(1, i);
    double v_abs = X_sigma_aug(2, i);
    double yaw = X_sigma_aug(3, i);
    double yawd = X_sigma_aug(4, i);
    double nu_a = X_sigma_aug(5, i);
    double nu_yawdd = X_sigma_aug(6, i);

    // predicted state values - Deterministic Part
    double px_p, py_p;

    if(fabs(yawd) > 0.001)
    {
      px_p = px + v_abs/yawd*(sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = py + v_abs/yawd*(-cos(yaw + yawd*delta_t) + cos(yaw));
    }
    else
    {
      px_p = px + v_abs* cos(yaw)*delta_t;
      py_p = py + v_abs* sin(yaw)*delta_t;
    }

    double v_p = v_abs;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // Add Stochastic Noise to the state
    px_p = px_p + 0.5*delta_t*delta_t*nu_a*cos(yaw);
    py_p = py_p + 0.5*delta_t*delta_t*nu_a*sin(yaw);
    v_p = v_p + nu_a*delta_t;
    yaw_p = yaw_p + 0.5*delta_t*delta_t*nu_yawdd;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //Predicted Sigma Points
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;

  }

  /*************************************************************
  *      Prediction mean and Covariance from Sigma Points     
  *************************************************************/
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for(int i =1; i<2*n_aug_+1; i++){
    double weight_tmp = 0.5/(lambda_+n_aug_);
    weights_(i) = weight_tmp;
  }  

  // State Mean
  x_.fill(0.0); 
  for(int i = 0; i < 2*n_aug_+1; i++){
    x_ = x_ + weights_(i)*Xsig_pred_.col(i);
  } 

  // State Covariance
  P_.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++){
    VectorXd xdiff = Xsig_pred_.col(i) - x_;    

    while (xdiff(3)> M_PI) xdiff(3) -= 2.*M_PI;
    while (xdiff(3)< -M_PI) xdiff(3) += 2.*M_PI;

    P_ = P_ + weights_(i) * xdiff*xdiff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  /*************************************************************
  *       Predicted Measurement Sigma Points     
  *************************************************************/
  VectorXd z = meas_package.raw_measurements_;

  // number of measured states using lidar px & py
  int n_s = 2; 

  MatrixXd Z_sigma = MatrixXd(n_s, 2 * n_aug_ + 1);

  for(int i =0; i < 2*n_aug_+1; i++){
    float px = Xsig_pred_(0,i);
    float py = Xsig_pred_(1,i);

    Z_sigma(0,i) = px;
    Z_sigma(1,i) = py;
  }

  /*************************************************************
  *        Predicted Measurements - Mean     
  *************************************************************/
  VectorXd z_pred = VectorXd(n_s);
  z_pred.fill(0.0); 
  for(int i = 0; i < 2*n_aug_+1; i++){
    z_pred = z_pred + weights_(i)*Z_sigma.col(i);
  } 

  /*************************************************************
  *       Predicted Measurement Covariance  - Mean    
  *************************************************************/  
  MatrixXd S = MatrixXd(n_s, n_s);

  S.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++){
    VectorXd zdiff = Z_sigma.col(i) - z_pred;

    S = S + weights_(i)*zdiff*zdiff.transpose();
  }

  MatrixXd R_laser_ = MatrixXd(n_s, n_s);

  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

  S = S + R_laser_;

  /*************************************************************
  *       UPDATE STEP: Calculate Cross Correlation Matrix    
  *************************************************************/
  MatrixXd T = MatrixXd(n_x_, n_s);
  T.fill(0.0);

  for(int i = 0; i < 2*n_aug_+1; i++){

    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    VectorXd zdiff = Z_sigma.col(i) - z_pred; 

    T = T + weights_(i)*xdiff*zdiff.transpose();
  }

  /*************************************************************
  *       UPDATE STEP: Kalman Gain    
  *************************************************************/

  MatrixXd K = T*S.inverse();

  /*************************************************************
  *       UPDATE STEP: State & Covariance Update  
  *************************************************************/
  // redidual 
  VectorXd y = z - z_pred;

  x_ = x_ + K * y;
  P_ = P_ - K*S*K.transpose();

  /*************************************************************
  *               NIS   
  *************************************************************/
  NIS_laser_ = y.transpose()* S.inverse() * y;

  NIS_ = NIS_laser_;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  /*************************************************************
  *       Predicted Measurement Sigma Points     
  *************************************************************/
  VectorXd z = meas_package.raw_measurements_;

  // number of measured states using Radar rho, phi, rho_dot
  int n_s = 3; 

  MatrixXd Z_sigma = MatrixXd(n_s, 2 * n_aug_ + 1);

  for(int i =0; i < 2*n_aug_+1; i++)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = v * cos(yaw);
    double v2 = v* sin(yaw);

    double ro  = sqrt(px*px + py*py);

    if (fabs(ro) < 0.001){
      ro = 0.001;
    }

    double phi = atan2(py, px);
    double ro_dot = (px*v1 + py*v2)/ro;

    Z_sigma(0,i) = ro;
    Z_sigma(1,i) = phi;
    Z_sigma(2,i) = ro_dot;
  }

  /*************************************************************
  *        Predicted Measurements - Mean     
  *************************************************************/
  VectorXd z_pred = VectorXd(n_s);
  z_pred.fill(0.0); 
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    z_pred = z_pred + weights_(i)*Z_sigma.col(i);
  } 

  /*************************************************************
  *       Predicted Measurement Covariance  - Mean    
  *************************************************************/  
  MatrixXd S = MatrixXd(n_s, n_s);

  S.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    VectorXd zdiff = Z_sigma.col(i) - z_pred;

    //normalization
    while(zdiff(1) > M_PI) zdiff(1) -= 2.*M_PI;
    while(zdiff(1) < -M_PI) zdiff(1) += 2.*M_PI;


    S = S + weights_(i)*zdiff*zdiff.transpose();
  }

  MatrixXd R_radar_ = MatrixXd(n_s, n_s);

  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0, std_radrd_*std_radrd_;

  S = S + R_radar_;

  /*************************************************************
  *       UPDATE STEP: Calculate Cross Correlation Matrix    
  *************************************************************/
  MatrixXd T = MatrixXd(n_x_, n_s);
  T.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {

    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    while(xdiff(1) > M_PI) xdiff(3) -= 2.*M_PI;
    while(xdiff(1) < -M_PI) xdiff(3) += 2.*M_PI;

    VectorXd zdiff = Z_sigma.col(i) - z_pred; 
    while(zdiff(1) > M_PI) zdiff(1) -= 2.*M_PI;
    while(zdiff(1) < -M_PI) zdiff(1) += 2.*M_PI;

    T = T + weights_(i)*xdiff*zdiff.transpose();
  }

  /*************************************************************
  *       UPDATE STEP: Kalman Gain    
  *************************************************************/

  MatrixXd K = T*S.inverse();

  /*************************************************************
  *       UPDATE STEP: State & Covariance Update  
  *************************************************************/
  // redidual 
  VectorXd y = z - z_pred;

  //angle normalization
  while(y(1) > M_PI) y(1) -= 2.*M_PI;
  while(y(1) < -M_PI) y(1) += 2.*M_PI;

  x_ = x_ + K * y;
  P_ = P_ - K*S*K.transpose();

  /*************************************************************
  *               NIS   
  *************************************************************/
  NIS_radar_ = y.transpose()* S.inverse() * y;

  NIS_ = NIS_radar_;

  
}

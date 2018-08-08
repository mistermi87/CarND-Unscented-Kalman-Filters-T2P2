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
UKF::UKF() {
  is_initialized_ = false;

  time_us_ = 0;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.4; //(~M_PI/8)

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

  // State dimension
  n_x_=x_.size();

  // Augmented state dimension
  n_aug_=n_x_+2;

  // Sigma point spreading parameter
  lambda_=3 - n_aug_;

  // Weights of sigma points
  weights_=VectorXd(2*n_aug_+1);
  //set weights
  weights_(0)=lambda_/(lambda_+n_aug_);
  for(int i=1;i<2*n_aug_+1;i++){
      weights_(i)=1.0/(2.0*(lambda_+n_aug_));
  };

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

   // agumented sigma points
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ +1 );

  //set covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  //initialize NIS counters
  total_las_meas_=0;
  NIS95_las_meas_=0;
  total_rad_meas_=0;
  NIS95_rad_meas_=0;

}
//Destructor
UKF::~UKF() {}


/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


  if (!is_initialized_) {
    // first measurement
    //cout << "EKF: " << endl;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        //cout << "1. Measurement Radar:"  << meas_package.raw_measurements_ << endl;

      ///*    Convert radar from polar to cartesian coordinates and initialize state.

      float ro = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      float ro_dot = meas_package.raw_measurements_(2);
      x_ <<ro*cos(theta), ro*sin(theta), sqrt(ro_dot), theta,0;

      double sigma_px= std_radr_*cos(std_radphi_);
      double sigma_py= std_radr_*sin(std_radphi_);

      double sigma2_px=sigma_px*sigma_px;
      double sigma2_py=sigma_py*sigma_py;
      double sigma2_pv= std_radrd_*std_radrd_;

      P_ <<  sigma2_px, 0, 0, 0, 0,
             0, sigma2_py, 0, 0, 0,
             0, 0, sigma2_pv, 0, 0,
             0, 0, 0, 1, 0,
             0, 0, 0, 0, 1;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        //cout << "1. Measurement Laser:"  << meas_package.raw_measurements_ << endl;
      x_ <<meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 4, 0.001, 0;

      double sigma2_px= std_laspx_*std_laspx_;
      double sigma2_py= std_laspy_*std_laspx_;

      P_ <<  sigma2_px, 0, 0, 0, 0,
             0, sigma2_py, 0, 0, 0,
             0, 0, 1, 0, 0,
             0, 0, 0, 1, 0,
             0, 0, 0, 0, 1;
    }


    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
      //std::cout << "is initialized, x_:" << std::endl<< x_<< std::endl;
      //std::cout << "is initialized, P:" << std::endl<< P_<< std::endl;
    return;
  }

  ///* Prediction:--------------------------------------

  //compute delta_t, reset previous timestamp for next delta calculation
  double delta_t= (meas_package.timestamp_-time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  ///* Update:--------------------------------------


  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
  }
    //std::cout << "Predicted state after Measurement" << std::endl;
    //std::cout << x_ << std::endl;
    //std::cout << "Predicted covariance matrix after Measurement" << std::endl;
    //std::cout << P_ << std::endl;

}







/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  ///* Generate Sigma Points:--------------------------------------
    //std::cout<<"===== Prediction ====="<<std::endl<<"--Generate Sigma-points--"<< std::endl;
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented mean state
  VectorXd x_aug =VectorXd(n_aug_);
  x_aug.head(n_x_)=x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_)=P_;
  P_aug.bottomRightCorner(2,2)<< std_a_*std_a_,0,
                                0, std_yawdd_*std_yawdd_;
    //std::cout<<"P_aug: "<<P_aug<< std::endl;

  //create square root matrix
  MatrixXd A=P_aug.llt().matrixL();
    //std::cout<<"A"<<A<< std::endl;

  //create augmented sigma points
  MatrixXd root_part;
  root_part=sqrt(lambda_+n_aug_)*A;

  for(int i=0;i<n_aug_*2+1;i++) {
   Xsig_aug_.col(i)=x_aug;
  };

  Xsig_aug_.block(0,1,n_aug_,n_aug_)+=root_part;
  Xsig_aug_.block(0,n_aug_+1, n_aug_,n_aug_)-=root_part;

    //std::cout << "Xsig_aug = " << Xsig_aug_ << std::endl;
    //std::cout << "Delta_t = " << delta_t << std::endl;



  ///*Predict Sigma Points ----------------------
   //std::cout<<"--Predict Sigma-points--"<< std::endl;
  for(int i=0; i<2 * n_aug_ + 1;i++){
      double px = Xsig_aug_(0,i);
      double py = Xsig_aug_(1,i);
      double v = Xsig_aug_(2,i);
      double yaw = Xsig_aug_(3,i);
      double yawd = Xsig_aug_(4,i);
      double nu_a = Xsig_aug_(5,i);
      double nu_yawdd = Xsig_aug_(6,i);

      double px_p, py_p;
      if (fabs(yawd) > 0.001) {
          px_p = px + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
          py_p = py + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
      }
      else {
          px_p = px + v*delta_t*cos(yaw);
          py_p = py + v*delta_t*sin(yaw);
      }
      px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
      py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
      //std::cout << "px_p, py_p" << std::endl;

      double v_p=v+delta_t*nu_a;

      double yaw_p= yaw+yawd*delta_t+0.5*delta_t*delta_t*nu_yawdd;

      double yawd_p=yawd+delta_t*nu_yawdd;

      // avoid zero division
      if (fabs(px_p) < 0.0001 && fabs(py_p) < 0.0001) px_p = 0.0001;

      Xsig_pred_.col(i)<<px_p, py_p, v_p, yaw_p, yawd_p;
        //std::cout << "Working on col" << i<< std::endl;

  };

    //std::cout << "Xsig_pred = " << Xsig_pred_ << std::endl;


  ///*predict state mean-----------------------------------
  for(int i=0;i<n_x_;i++){
     x_(i)=Xsig_pred_.row(i)*weights_;
  };
     //std::cout << "Predicted state" << std::endl;
     //std::cout << x_ << std::endl;

  ///*predict state covariance matrix------------------------------
  P_.fill(0.0);
  for(int i=0;i<2 * n_aug_ + 1;i++){
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  };
    //std::cout << "Predicted covariance matrix" << std::endl;
    //std::cout << P_ << std::endl;
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    //std::cout << "===== Update Lidar =====" << std::endl;
    int n_z = 2;
  //create matrix for sigma points in measurement space

  //MatrixXd Zsig = Xsig_pred_.topRows<2>();
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, n_aug_*2+1);

  //mean predicted measurement
  VectorXd z_pred = Zsig*weights_;

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  //calculate innovation covariance matrix S
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<   std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;


  for(int i=0;i<2 * n_aug_ + 1;i++){
    VectorXd diff=Zsig.col(i)-z_pred;

    S=S+ weights_(i)*diff*diff.transpose();
  }
  S=S+R;
    //std::cout << "sigma points in measurement space Zsig" << std::endl << Zsig << std::endl;
    //std::cout << "z_pred" << std::endl << z_pred << std::endl;
    //std::cout << "measurement covariance matrix S" << std::endl << S << std::endl;

  ///*Update----------------------------------------------------------

  //Get data of current measurement
  VectorXd z=meas_package.raw_measurements_;
  //std::cout << "Raw measuremnts z" << std::endl << z << std::endl;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for(int i=0;i<2*n_aug_+1;i++){
      VectorXd diff1 = Xsig_pred_.col(i) - x_;
      while(diff1(3)>M_PI) diff1(3)-=2*M_PI;
      while(diff1(3)<-M_PI) diff1(3)+=2*M_PI;

      VectorXd diff2 = Zsig.col(i) - z_pred;

      Tc=Tc+weights_(i)*diff1*diff2.transpose();
  };
    //std::cout << "cross correlation matrix Tc" << std::endl;
    //std::cout << Tc << std::endl;

  //calculate Kalman gain K;
  MatrixXd K = Tc* S.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff = z-z_pred;
  x_=x_+ K * z_diff;

  P_=P_-K*S*K.transpose();

  ///* NIS Calculation
/*   float eps = z_diff.transpose()*S.inverse()*z_diff;
 *   total_las_meas_++;
 *   if(eps>5.991){
 *     NIS95_las_meas_++;
 *   }
 *    std::cout <<"NIS Lidar "<<  eps << "\t || " << NIS95_las_meas_ << " of " << total_las_meas_ << " measurements are above 95% for chi2 distr." << std::endl;
 */

}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    //std::cout << "===== Update Radar =====" << std::endl;

  int n_z = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  //transform sigma points into measurement space
  for(int i=0;i<2 * n_aug_ + 1;i++){
      Zsig(0,i)=sqrt((Xsig_pred_(0,i))*(Xsig_pred_(0,i))+(Xsig_pred_(1,i))*(Xsig_pred_(1,i)));
      Zsig(1,i)=atan2((Xsig_pred_(1,i)),(Xsig_pred_(0,i)));
      Zsig(2,i)=(Xsig_pred_(0,i)*cos(Xsig_pred_(3,i))*Xsig_pred_(2,i)+Xsig_pred_(1,i)*sin(Xsig_pred_(3,i))*Xsig_pred_(2,i))/Zsig(0,i);
  }
    //std::cout << "sigma points in measurement space Zsig" << std::endl << Zsig << std::endl;
  //calculate mean predicted measurement
  z_pred=Zsig*weights_;
    //std::cout << "z_pred" << std::endl << z_pred << std::endl;
  //calculate innovation covariance matrix S
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<   std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_,0,
        0,0,std_radrd_*std_radrd_;

  for(int i=0;i<2 * n_aug_ + 1;i++){
    VectorXd diff=Zsig.col(i)-z_pred;

    while (diff(1)> M_PI) diff(1)-=2.*M_PI;
    while (diff(1)< -M_PI) diff(1)+=2.*M_PI;

    S=S+ weights_(i)*diff*diff.transpose();
  };
  S=S+R;
    //std::cout << "measurement covariance matrix S" << std::endl << S << std::endl;

  ///*Update----------------------------------------------------------

  //Get data of current measurement
  VectorXd z=meas_package.raw_measurements_;
    //std::cout << "Actual measurement z" << std::endl << z << std::endl;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  //calculate cross correlation matrix
  for(int i=0;i<2*n_aug_+1;i++){
      VectorXd diff1 = Xsig_pred_.col(i) - x_;
      while(diff1(3)>M_PI) diff1(3)-=2*M_PI;
      while(diff1(3)<-M_PI) diff1(3)+=2*M_PI;

      VectorXd diff2 = Zsig.col(i) - z_pred;
      while(diff2(1)>M_PI) diff2(1)-=2*M_PI;
      while(diff2(1)<-M_PI) diff2(1)+=2*M_PI;

      Tc=Tc+weights_(i)*diff1*diff2.transpose();
  };
     //std::cout << "cross correlation matrix Tc" << std::endl << Tc << std::endl;
  //calculate Kalman gain K;
  MatrixXd K = Tc* S.inverse();
  //update state mean and covariance matrix
  //std::cout << "Kalman gain" << std::endl << K << std::endl;
  VectorXd z_diff = z-z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    //std::cout << "z_diff" << std::endl << z_diff << std::endl;

  x_=x_+ K * z_diff;


  P_=P_-K*S*K.transpose();

  ///* NIS Calculation
/*   float eps = z_diff.transpose()*S.inverse()*z_diff;
 *   total_rad_meas_++;
 *   if(eps>7.815){
 *     NIS95_rad_meas_++;
 *   }
 *   std::cout <<"NIS Radar "<<  eps << "\t || " << NIS95_rad_meas_ << " of " << total_rad_meas_ << " measurements are above 95% for chi2 distr." << std::endl;
 */

}

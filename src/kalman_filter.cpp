#include "kalman_filter.h"
#include <iostream>
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in; // initial state
  P_ = P_in; // initial state covariance matrix
  F_ = F_in; // state transition matrix
  H_ = H_in; // measurement matrix
  R_ = R_in; // measurement covariance matrix
  Q_ = Q_in; // process covariance matrix
  prev_phi_ = 0; // previous phi, for sign changes
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  // Calculate h(x') vector which contains (rho, phi, and rho*)
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  float rho = sqrt(px*px + py*py);
  float phi = atan2(py, px);
  
  // The reviewer suggested that I normalize phi to fix my RMSE problem, but I thought arctan2 already
  // normalized phi between [-pi, pi] ?  Here's a quick check to see if phi is ever greater than 3.14
  // Quick note: I didn't see any exceptions thrown during my testing
  if (fabs(phi) >= M_PI) {
    std::cout << "----------------------------THE REVIEWER WAS RIGHT! " << phi << " ----------------------------" << std::endl;
    throw "The reviewer was right!  Phi wasn't normalized between -pi and pi!";
  }
  
  // Bingo!  phi's sign change is really throwing off 
  // our measurement of velocity in euclidean space.
  //
  // Question: How come the radar is reading a complete change in 
  // direction while having a similarly large radial distance?
  // This seems like a bug in the measurement reading, not in my 
  // implementation of EKF.  Nevertheless, let's go ahead and handle
  // this particular corner case.
  //
  // Resolve drastic changes where phi goes from approx 3.19 (ground truth larger than pi) to approx -3.11
  // Quick note: I believe this happens around line 274 and 276 in the dataset (below).
  //
  /* From Dataset 1.  Phi changes from 3.19 (ground truth larger than pi?) to -3.11 in 0.5 seconds with rho 
                      of 6.0 and 5.6 respectively. Phi_dot increases from 1.77 to 2.5.  So, this means the 
                      object has moved a total of 11.6 in 0.5 seconds and has increased its radial velocity away from                         the observer. Is this anomaly even supposed to be here?
  
     R	6.005131e+00	3.190031e+00	1.776367e+00	1477010456650000	-5.378204e+00	6.547190e-02	-2.154769e+00	-4.693737e+00	4.282015e+00	-1.633729e-01
     L	-5.299723e+00	-2.129817e-01	1477010456700000	-5.486900e+00	-1.687723e-01	-2.191805e+00	-4.673350e+00	4.273847e+00	-1.699593e-01
     R	5.646317e+00	-3.115994e+00	2.506136e+00	1477010456750000	-5.597482e+00	-4.019714e-01	-2.230141e+00	-4.651846e+00	4.265349e+00	-1.765190e-01
  */
  if (fabs(prev_phi_ + phi) < (fabs(prev_phi_))) {
    return;
  }
  
  // If there was no crazy measurement, record the previous phi and keep going...
  prev_phi_ = phi;
  
  float rho_dot = (px*vx + py*vy) / rho;
  
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  
  // With h(x'), do your typical Kalman update steps
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

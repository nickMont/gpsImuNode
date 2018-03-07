#include "filterTW.h"
#include <Eigen/LU>  // For matrix inverse

namespace gpsimu_odom {

void KalmanTW::initialize(const State_t &state,
                              const StateCov_t &initial_cov,
                              const ProcessCov_t &process_noise,
                              const MeasurementCov_t &meas_noise) {
  x = state;
  P = initial_cov;
  Q = process_noise;
  R = meas_noise;
}

void KalmanTW::processUpdate(double dt, Eigen::Vector3d uT) {
  //A, 7x7
  StateCov_t A = StateCov_t::Identity();
  A(0,3) = dt; A(1,4)=dt; A(2,5)=dt;
  //B, 7x3
  Bmat_t B = Bmat_t::Zero();
  B.topRows(3) = dt*dt*0.5*Eigen::Matrix3d::Identity();
  B.middleRows(3,3) = Eigen::Matrix3d::Identity()*dt;
  //B*u*tw% is the effect of tw% on x(0-5) so should be the rightmost column of A(7x7)
  State_t A_rightCol = B*uT;
  A.rightCols(1) = A_rightCol;
  //noise matrix
  Gamma_t gamma = Gamma_t::Zero();
  gamma.topLeftCorner(3,3)=dt*dt*0.5*Eigen::Matrix3d::Identity();
  gamma.block(3,0,3,3)=dt*Eigen::Matrix3d::Identity();
  gamma(6,3)=1;

  //propagation
  x = A * x - B * Eigen::Vector3d(0,0,-9.81);
  P = A * P * A.transpose() + gamma*Q*gamma.transpose();
}

void KalmanTW::measurementUpdate(const Measurement_t &meas, double dt) {
  Eigen::Matrix<double, 3, 7> H;
  H.setZero();
  H.leftCols(3) = Eigen::Matrix3d::Identity();

  const Eigen::Matrix<double, 7, 3> K =
      P * H.transpose() * (H * P * H.transpose() + R).inverse();
  const Measurement_t inno = meas - H * x;
  x += K * inno;
  P = (StateCov_t::Identity() - K * H) * P;
}

}  // namespace gps_odom

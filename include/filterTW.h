#include <Eigen/Core>

namespace gpsimu_odom {

class KalmanTW {
 public:
  typedef Eigen::Matrix<double, 7, 1> State_t;
  typedef Eigen::Matrix<double, 7, 7> StateCov_t;
  typedef Eigen::Matrix<double, 4, 4> ProcessCov_t;
  typedef Eigen::Matrix<double, 3, 1> Measurement_t;
  typedef Eigen::Matrix<double, 3, 3> MeasurementCov_t;
  typedef Eigen::Matrix<double, 7, 3> Bmat_t;
  typedef Eigen::Matrix<double, 7, 4> Gamma_t;

  KalmanTW() {}
  KalmanTW(const State_t &state, const StateCov_t &initial_cov,
               const ProcessCov_t &process_noise,
               const MeasurementCov_t &meas_noise)
      : x(state), P(initial_cov), Q(process_noise), R(meas_noise) {}

  void initialize(const State_t &state, const StateCov_t &initial_cov,
                  const ProcessCov_t &process_noise,
                  const MeasurementCov_t &meas_noise);

  void processUpdate(double dt, Eigen::Vector3d uu);
  void measurementUpdate(const Measurement_t &meas, double dt);
  void setState(const State_t &new_state) {x = new_state;}
  void setEstimateCovariance(const StateCov_t &new_covariance) { P = new_covariance; }
  void setProcessNoise(const ProcessCov_t &process_noise) { Q = process_noise; }
  void setMeasurementNoise(const MeasurementCov_t &meas_noise) {
    R = meas_noise;
  }

  const State_t &getState() const { return x; }
  const StateCov_t &getCovariance() const { return P; }

 private:
  State_t x;
  StateCov_t P;
  ProcessCov_t Q;
  MeasurementCov_t R;
};

}  // namespace gps_odom

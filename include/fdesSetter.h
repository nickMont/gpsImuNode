#ifndef GPS_ODOM_FILTER_H_
#define GPS_ODOM_FILTER_H_

#include <Eigen/Core>

namespace gps_odom {

class fdesControl {
 public:
  typedef Eigen::Matrix<double, 3, 1> Fdes_t;

  const Fdes_t &getFdes() const { return uF; }
  void setFdes(utest) const{uF=utest;}
  void initialize(u0) const{uF=u0;}

 private:
  Fdes_t uF;
};

}  // namespace gps_odom

#endif  // GPS_ODOM_FILTER_H_

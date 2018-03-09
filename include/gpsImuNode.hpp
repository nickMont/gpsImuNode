#include <ros/ros.h>
#include <tf2_ros/transform_broadcaster.h>
#include <std_msgs/Float64.h>
#include <geometry_msgs/Pose.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>
#include <Eigen/Geometry>
#include <gbx_ros_bridge_msgs/SingleBaselineRTK.h>
#include <gbx_ros_bridge_msgs/Attitude2D.h>
#include <gbx_ros_bridge_msgs/Imu.h>
#include <gbx_ros_bridge_msgs/ImuConfig.h>
#include <gbx_ros_bridge_msgs/NavigationSolution.h>
#include <stdint.h>

#include "filter.h"
#include "filterTW.h"

namespace gpsimu_odom
{
class gpsImuNode
{
 public:
  gpsImuNode(ros::NodeHandle &nh);

  //void gpsCallback(const geometry_msgs::PoseStamped::ConstPtr &msg);
  void singleBaselineRTKCallback(const gbx_ros_bridge_msgs::SingleBaselineRTK::ConstPtr &msg);
  void attitude2DCallback(const gbx_ros_bridge_msgs::Attitude2D::ConstPtr &msg);
  void imuConfigCallback(const gbx_ros_bridge_msgs::ImuConfig::ConstPtr &msg);
  void imuDataCallback(const gbx_ros_bridge_msgs::Imu::ConstPtr &msg);
  void tOffCallback(const gbx_ros_bridge_msgs::NavigationSolution::ConstPtr &msg);
  void publishOdomAndMocap();
  Eigen::Matrix<double,15,15> getFmatrixCF(const double dt, const Eigen::Vector3d fB,
    const Eigen::Vector3d omegaB, const Eigen::Matrix3d RR);
  Eigen::Matrix<double,15,6> getGammakmatrixCF(const double dt, const Eigen::Matrix3d RR);
  Eigen::Matrix<double,3,15> getHkmatrixOneAntennaCF(const Eigen::Vector3d Lab, const Eigen::Matrix3d RR);
  Eigen::Matrix<double,6,15> getHkmatrixTwoAntennaCF(const Eigen::Vector3d Limu,
    const Eigen::Vector3d Ls2p, const Eigen::Matrix3d RR);
  void kfCFPropagate(const double dt0, const Eigen::Matrix<double,15,15> P0,
    const Eigen::Matrix<double,15,1> x0, const Eigen::Matrix<double,15,15> F0,
    const Eigen::Matrix3d RR, const Eigen::Matrix<double,6,6> Qk,
    Eigen::Matrix<double,15,15> Pbar, Eigen::Matrix<double,15,1> xBar);
  void kfCFMeasure2(const Eigen::Matrix<double,15,1> xBar, const Eigen::Vector3d drI, const Eigen::Vector3d drS2P,
    const Eigen::Matrix3d RR, const Eigen::Matrix<double,15,15> Pbar,
    const Eigen::Vector3d Limu, const Eigen::Vector3d Ls2p, const Eigen::Matrix<double,6,6> Rbar,
    Eigen::Matrix<double,15,15> Pkp1, Eigen::Matrix<double,15,1> xkp1);
  void runCF(const double dt0);
  Eigen::Matrix3d updateRBIfromGamma(const Eigen::Matrix3d R0, const Eigen::Vector3d gamma);
  Eigen::Matrix3d hatmat(const Eigen::Vector3d v1);
  Eigen::Matrix3d rotMatFromEuler(Eigen::Vector3d ee);
  Eigen::Matrix3d rotMatFromQuat(Eigen::Quaterniond qq);
  Eigen::Matrix3d rotMatFromWahba(const Eigen::VectorXd weights, const::Eigen::MatrixXd vI,
    const::Eigen::MatrixXd vB);
  Eigen::Quaterniond rotmat2quat(const Eigen::Matrix3d RR);
  Eigen::Vector3d unit3(const Eigen::Vector3d v1);

 private:
  void PublishTransform(const geometry_msgs::Pose &pose,
                        const std_msgs::Header &header,
                        const std::string &child_frame_id);

  ros::Publisher odom_pub_;
  ros::Publisher localOdom_pub_;
  ros::Publisher mocap_pub_;
  ros::Publisher internalPosePub_;
  std::string child_frame_id_;
  tf2_ros::TransformBroadcaster tf_broadcaster_;
  Eigen::Vector3d baseECEF_vector, baseENU_vector, WRW0_ecef, arenaRefCenter,
      internal_rI, internal_rC, internal_rImu, rRefImu, n_err, imuAccelMeas,
      imuAttRateMeas, l_imu, l_s2p;

  Eigen::Matrix3d Recef2enu, Rwrw, R_G2wrw, RBI;
  Eigen::Matrix<double,21,3> rCtildeCalib, rBCalib;
  Eigen::Matrix<double,15,1> xState;
  Eigen::Matrix<double,15,15> Fimu, Pimu;
  Eigen::Matrix<double,6,6> Qimu, Rk;

  bool publish_tf_;
  ros::Subscriber gps_sub_, rtkSub_, a2dSub_, imuSub_, imuConfigSub_, tOffsetSub_;
  //geometry_msgs::PoseStamped::ConstPtr initPose_;
  geometry_msgs::PoseStamped initPose_;
  geometry_msgs::PoseStamped centerInENU;
  //geometry_msgs::PoseStamped initPose_;
  Eigen::Quaterniond internalQuat, quaternionSetpoint;
  int centerFlag, internalSeq, sec_in_week;
  double lastRTKtime, lastA2Dtime, minTestStat, max_accel, throttleSetpoint, throttleMax,
      imuConfigAccel, imuConfigAttRate, tMeasOffset, pi, tLastProcessed;
  bool validRTKtest, validA2Dtest, kfInit, hasAlreadyReceivedA2D, hasAlreadyReceivedRTK, hasRBI;

  int32_t trefWeek, trefSecOfWeek;
  float trefFracSecs;
  uint64_t sampleFreqNum, sampleFreqDen;

  uint64_t one, imuTimeTrunc;

};

} // gps_odom

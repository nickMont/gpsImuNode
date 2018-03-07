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
#include <gbx_ros_bridge_msgs/ObservablesMeasurementTime.h>
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
  void tOffCallback(const gbx_ros_bridge_msgs::ObservablesMeasurementTime::ConstPtr &msg);
  Eigen::Matrix<double,15,15> getFmatrixCF(const double dt, const Eigen::Vector3d fB,
    const Eigen::Vector3d omegaB, const Eigen::Matrix3d RBI);
  Eigen::Matrix<double,15,6> getGammakmatrixCF(const double dt, const Eigen::Matrix3d RBI);
  Eigen::Matrix<double,3,15> getHkmatrixCF(const Eigen::Vector3d Lab, const Eigen::Matrix3d RBI);
  Eigen::Matrix3d hatmat(const Eigen::Vector3d v1);
  Eigen::Matrix3d rotMatFromEuler(Eigen::Vector3d ee);
  Eigen::Matrix3d rotMatFromQuat(Eigen::Quaterniond qq);

 private:
  void PublishTransform(const geometry_msgs::Pose &pose,
                        const std_msgs::Header &header,
                        const std::string &child_frame_id);

  KalmanFilter kf_;
  KalmanTW kfTW_;
  //gnssimu gpsimu_;
  KalmanFilter::State_t xCurr;
  ros::Publisher odom_pub_;
  ros::Publisher localOdom_pub_;
  ros::Publisher mocap_pub_;
  ros::Publisher internalPosePub_; //publishes /Valkyrie/pose to itself
  std::string child_frame_id_;
  tf2_ros::TransformBroadcaster tf_broadcaster_;
  Eigen::Vector3d baseECEF_vector, baseENU_vector, WRW0_ecef, arenaRefCenter,
      internalPose, n_err, imuAccelMeas, imuAttRateMeas;
  Eigen::Matrix3d Recef2enu;
  Eigen::Matrix<double,15,1> xState;
  Eigen::Matrix<double,15,15> Fimu, Pimu;
  bool publish_tf_;
  ros::Subscriber gps_sub_, rtkSub_, a2dSub_, imuSub_, imuConfigSub_, tOffsetSub_;
  //geometry_msgs::PoseStamped::ConstPtr initPose_;
  geometry_msgs::PoseStamped initPose_;
  geometry_msgs::PoseStamped centerInENU;
  //geometry_msgs::PoseStamped initPose_;
  Eigen::Quaterniond internalQuat, quaternionSetpoint;
  int centerFlag, internalSeq, sec_in_week;
  double lastRTKtime, lastA2Dtime, minTestStat, max_accel, throttleSetpoint, throttleMax,
      imuConfigAccel, imuConfigAttRate, tMeasOffset, pi;
  bool validRTKtest, validA2Dtest, kfInit, hasAlreadyReceivedA2D, hasAlreadyReceivedRTK;

  int32_t trefWeek, trefSecOfWeek;
  float trefFracSecs;
  uint64_t sampleFreqNum, sampleFreqDen;

  uint64_t one, imuTimeTrunc;

};

} // gps_odom

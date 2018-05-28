#include <ros/ros.h>
#include <tf2_ros/transform_broadcaster.h>
#include <std_msgs/Float64.h>
#include <geometry_msgs/Pose.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <gbx_ros_bridge_msgs/SingleBaselineRTK.h>
#include <gbx_ros_bridge_msgs/Attitude2D.h>
#include <gbx_ros_bridge_msgs/Imu.h>
#include <gbx_ros_bridge_msgs/ImuConfig.h>
#include <gbx_ros_bridge_msgs/NavigationSolution.h>
#include <gbx_ros_bridge_msgs/ObservablesMeasurementTime.h>
//#include <stdio.h>
//#include <stdint.h>
#include <cmath>

#include "filter.h"
#include "filterTW.h"
//#include "classes.h" //Included in filterImu.hpp
#include "filterImu.hpp"

namespace gpsimu_odom
{
class estimationNode
{
 public:
    estimationNode(ros::NodeHandle &nh);

    //ROS stuff
    void singleBaselineRTKCallback(const gbx_ros_bridge_msgs::SingleBaselineRTK::ConstPtr &msg);
    void attitude2DCallback(const gbx_ros_bridge_msgs::Attitude2D::ConstPtr &msg);
    void imuConfigCallback(const gbx_ros_bridge_msgs::ImuConfig::ConstPtr &msg);
    void imuDataCallback(const gbx_ros_bridge_msgs::Imu::ConstPtr &msg);
    void navsolCallback(const gbx_ros_bridge_msgs::NavigationSolution::ConstPtr &msg);
    void tOffsetCallback(const gbx_ros_bridge_msgs::ObservablesMeasurementTime::ConstPtr &msg);
    void publishOdomAndMocap();

    //Math helper functions
    Eigen::Matrix3d updateRBIfromGamma(const Eigen::Matrix3d R0, const Eigen::Vector3d gamma);
    Eigen::Matrix3d hatmat(const Eigen::Vector3d v1);
    Eigen::Matrix3d rotMatFromWahba(const Eigen::VectorXd weights, const::Eigen::MatrixXd &vI,
        const::Eigen::MatrixXd &vB);
    Eigen::Quaterniond rotmat2quat(const Eigen::Matrix3d RR);
    Eigen::Vector3d unit3(const Eigen::Vector3d v1);
    Eigen::Matrix3d orthonormalize(const Eigen::Matrix3d inmat);
    Eigen::Matrix3d rotMatFromQuat(const Eigen::Quaterniond qq);
    void saturateBiases(const double baMax, const double bgMax);
    double symmetricSaturationDouble(const double inval, const double maxval);
    Eigen::Matrix3d euler2dcm312(const Eigen::Vector3d ee);
    Eigen::Matrix3d euler2dcm321(Eigen::Vector3d ee);

 private:
    void PublishTransform(const geometry_msgs::Pose &pose,
                                                const std_msgs::Header &header,
                                                const std::string &child_frame_id);

    ros::Publisher localOdom_pub_, mocap_pub_;
    std::string child_frame_id_;
    gpsImu imuFilter_;
    imuMeas lastImuMeas_;

    tf2_ros::TransformBroadcaster tf_broadcaster_;

    Eigen::Vector3d zeroInECEF_, zeroInWRW_, rPrimaryMeas_, rS2PMeas_, rPrimaryMeas_mu, Lcg2p_, Ls2p_, Lcg2imu_;
    Eigen::Matrix3d Recef2enu_, Rwrw_, Recef2wrw_, RBI_, QgyroOutput_;
    Eigen::Matrix<double,21,3> rCtildeCalib_, rBCalib_;

    ros::Subscriber gps_sub_, rtkSub_, a2dSub_, imuSub_, imuConfigSub_, tOffsetSub_, navSub_;

    int internalSeq;
    double lastRTKtime_, lastA2Dtime_, minTestStat_, imuConfigAccel_, imuConfigAttRate_, 
        pi, tLastProcessed_, tnavsolFracSecs_, sec_in_week_, toffsetFracSecs_, dtRXinMeters_;
    bool validRTKtest_, validA2Dtest_, hasAlreadyReceivedA2D_, hasAlreadyReceivedRTK_, rbiIsInitialized_,
        isCalibrated_, publish_tf_;

    long long int tIndexConfig_;
    uint32_t toffsetWeek_, toffsetSecOfWeek_;
    uint64_t sampleFreqNum_, sampleFreqDen_;

    std::string updateType;

};

} // gps_odom

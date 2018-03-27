#include <Eigen/Geometry>
#include "gpsImuNode.hpp"
#include <string>
#include <iostream>

Eigen::Matrix3d ecef2enu_rotMatrix(Eigen::Vector3d ECEF){


    //----- Define WGS-84 Earth parameters
    const double aa = 6378137.00000;
    const double bb = 6356752.3142518;
    const double ee = 0.0818191908334158;
    const double ep = sqrt((aa*aa - bb*bb)/(bb*bb));
    const double ee2 = (aa*aa-bb*bb)/(aa*aa);

    //----- Convert to (phi,lambda,h) geodetic coordinates
    double x = ECEF(0);
    double y = ECEF(1);
    double z = ECEF(2);
    double lambda = atan2(y, x);
    double p = sqrt(x*x + y*y);
    double theta = atan2(z*aa, p*bb);
    /*double phi = atan2(z + ep*ep*bb*pow(sin(theta),3),
                    p - ee*ee*aa*pow(cos(theta),3));*/
    double phi=atan2(z,(1-ee2)*p);
    double N,h, phiprev;
    bool contvar=true;
    while(contvar)
    {
        phiprev=phi;
        N=aa/sqrt(1-ee2*sin(phi)*sin(phi));
        h=p/cos(phi)-N;
        phi=atan2(z,(1-ee2*N/(N+h))*p);
        if(abs(phiprev-phi)<1e-6)
        {
            contvar=false;
        }
    }

    //----- Form the rotation matrix
    Eigen::Matrix3d Renu_ecef = Eigen::Matrix3d::Zero();
    Renu_ecef(0,0) = -sin(lambda);
    Renu_ecef(0,1) = cos(lambda);
    Renu_ecef(0,2) = 0;
    Renu_ecef(1,0) = -sin(phi)*cos(lambda);
    Renu_ecef(1,1) = -sin(phi)*sin(lambda);
    Renu_ecef(1,2) = cos(phi);
    Renu_ecef(2,0) = cos(phi)*cos(lambda);
    Renu_ecef(2,1) = cos(phi)*sin(lambda);
    Renu_ecef(2,2) = sin(phi);

 return Renu_ecef;
}

Eigen::Vector3d ecef2enu(Eigen::Vector3d ECEF){
    Eigen::Matrix3d R = ecef2enu_rotMatrix(ECEF);
    Eigen::Vector3d ENU = R*ECEF;
    /*ROS_INFO("R row 1: %f %f %f", R(0,0), R(0,1), R(0,2));
    ROS_INFO("handmat: %f %f %f", R(0,0)*ECEF(0),  R(0,1)*ECEF(1),  R(0,2)*ECEF(2));
    ROS_INFO("ENUy: %f  ECEFy: %f", ENU(1), ECEF(1)); */
    return ENU;
}

namespace gpsimu_odom
{
gpsImuNode::gpsImuNode(ros::NodeHandle &nh)
{

  //Get data about node and topic to listen
  std::string quadPoseTopic, quadName, rtktopic, a2dtopic, posePubTopic, nodeNamespace;
  double tmax;
  quadName = ros::this_node::getName();
//  Eigen::Vector3d enuInput;
  nodeNamespace = ros::this_node::getNamespace();    
  ros::param::get(quadName + "/quadPoseTopic", quadPoseTopic);
  ros::param::get(quadName + "/arenaCenterX", baseECEF_vector(0));
  ros::param::get(quadName + "/arenaCenterY", baseECEF_vector(1));
  ros::param::get(quadName + "/arenaCenterZ", baseECEF_vector(2));
  ros::param::get(quadName + "/arenaCenterX_ENU", n_err(0));
  ros::param::get(quadName + "/arenaCenterY_ENU", n_err(1));
  ros::param::get(quadName + "/arenaCenterZ_ENU", n_err(2));
  ros::param::get(quadName + "/rtktopic", rtktopic);
  ros::param::get(quadName + "/a2dtopic", a2dtopic);
  ros::param::get(quadName + "/posePubTopic", posePubTopic);
  ros::param::get(quadName + "/minimumTestStat",minTestStat);
  ros::param::get(quadName + "/maxThrust",tmax);

  Eigen::Matrix<double,6,6> temp;
  Qimu=1e-1*Eigen::Matrix<double,6,6>::Identity();
  Rk=1e-2*Eigen::Matrix<double,6,6>::Identity();

  one = 1UL;

  //Get additional parameters for the kalkman filter
  nh.param(quadName + "/max_accel", max_accel, 2.0);
  nh.param(quadName + "/publish_tf", publish_tf_, true);
  nh.param<std::string>(quadName + "/child_frame_id", child_frame_id_, "base_link");
  if(publish_tf_ && child_frame_id_.empty())
    throw std::runtime_error("gpsimu_odom: child_frame_id required for publishing tf");

  //should be a const but catkin doesn't like scoping it
  pi = std::atan(1.0)*4;

  Recef2enu=ecef2enu_rotMatrix(baseECEF_vector);

  //Account for WRW rotation wrt ENU
  double thetaWRW;
  thetaWRW = 6.2*pi/180; //angle of rooftop coordinate system WRT ENU
  Rwrw << cos(thetaWRW), -1*sin(thetaWRW), 0,
          sin(thetaWRW), cos(thetaWRW), 0,
          0, 0, 1;
  Pimu=Eigen::MatrixXd::Identity(15,15);
  R_G2wrw=Rwrw*Recef2enu;

  //THIS SHOULD BE INITIALIZED IN GPSCALLBACKS
  rRefImu<<0,0,0; //location of rI_0 for imu

  lastRTKtime=0;
  lastA2Dtime=0;
  //internalQuat.resize(4);
  internalSeq=0;

  //Secondary to primary vector in body frame
  l_s2p<<0.195,0,0;
  l_imu<<-0.195,0,-0.10;

  hasRBI = false;

  // Initialize publishers and subscribers
  //odom_pub_ = nh.advertise<nav_msgs::Odometry>("odom", 10); //MUST have a node namespace, ns="quadName", in launchfile
  localOdom_pub_ = nh.advertise<nav_msgs::Odometry>("local_odom_INS", 10);
  mocap_pub_ = nh.advertise<geometry_msgs::PoseStamped>("mavros/mocap/pose", 10);
/*  gps_sub_ = nh.subscribe(quadPoseTopic, 10, &gpsImuNode::gpsCallback,
                            this, ros::TransportHints().tcpNoDelay());*/
//  internalPosePub_ = nh.advertise<geometry_msgs::PoseStamped>(posePubTopic,10);
  rtkSub_ = nh.subscribe("SingleBaselineRTK",10,&gpsImuNode::singleBaselineRTKCallback,
                            this, ros::TransportHints().tcpNoDelay());
  a2dSub_ = nh.subscribe("Attitude2D",10,&gpsImuNode::attitude2DCallback,
                            this, ros::TransportHints().tcpNoDelay());
  imuSub_ = nh.subscribe("IMU",10, &gpsImuNode::imuDataCallback,
                            this, ros::TransportHints().tcpNoDelay());
  imuConfigSub_ = nh.subscribe("IMUConfig",10, &gpsImuNode::imuConfigCallback,
                            this, ros::TransportHints().tcpNoDelay());
  tOffsetSub_ = nh.subscribe("NavigationSolution",10,&gpsImuNode::tOffCallback,
                            this, ros::TransportHints().tcpNoDelay());
  ROS_INFO("Waiting for IMU config data, this may take a moment...");
  gbx_ros_bridge_msgs::ImuConfig::ConstPtr imuConfigMsg = ros::topic::waitForMessage<gbx_ros_bridge_msgs::ImuConfig>("IMUConfig");
  imuConfigAccel = imuConfigMsg->lsbToMetersPerSecSq;
  imuConfigAttRate = imuConfigMsg->lsbToRadPerSec;
  sampleFreqNum = imuConfigMsg->sampleFreqNumerator;
  sampleFreqDen = imuConfigMsg->sampleFreqDenominator;  
  tIndexKconfig = imuConfigMsg->tIndexk;
  ROS_INFO("IMU configuration recorded, finishing startup.");
  ROS_INFO("Startup complete");

  //Get initial pose
  //initPose_ = ros::topic::waitForMessage<geometry_msgs::PoseStamped>(quadPoseTopic);
  //geometry_msgs::PoseStamped initPose_;
}


void gpsImuNode::tOffCallback(const gbx_ros_bridge_msgs::NavigationSolution::ConstPtr &msg)
{
  //static int SEC_PER_WEEK = 604800;
  //tMeasOffset=msg->tOffset.week*SEC_PER_WEEK+msg->tOffset.secondsOfWeek+msg->tOffset.fractionOfSecond;
  //time offset from converting RRT to ORT
  trefWeek = msg->tSolution.week;
  trefFracSecs = msg->tSolution.fractionOfSecond;
  trefSecOfWeek = msg->tSolution.secondsOfWeek;
}


void gpsImuNode::imuConfigCallback(const gbx_ros_bridge_msgs::ImuConfig::ConstPtr &msg)
{
  std::cout << "Config message received!" << std::endl;
  imuConfigAccel = msg->lsbToMetersPerSecSq; //scaling to m/s2 from "non-engineering units"
  imuConfigAttRate = msg->lsbToRadPerSec; //scaling to rad/s from "non-engineering units"
  //imuSampleFreq = msg->sampleFreqNumerator/msg->sampleFreqDenominator/36/3600;  //samples per second
  
  sampleFreqNum = msg->sampleFreqNumerator;
  sampleFreqDen = msg->sampleFreqDenominator;
  tIndexKconfig = msg->tIndexk;
}


void gpsImuNode::imuDataCallback(const gbx_ros_bridge_msgs::Imu::ConstPtr &msg)
{
  //Initialization variables
  static bool isCalibrated(false);
  static int counter(0);
  static Eigen::Matrix<double,100,3> imuAStore, imuGStore;
  static float tLastImu(0);
  static Eigen::Vector3d ba0(0,0,0);
  static Eigen::Vector3d bg0(0,0,0);
  //gps variables
  // static const int SEC_PER_WEEK(604800);
  static const int SF_TL = 24;
  static const int32_t SF_T = 0x1 << SF_TL;
  static const int32_t SF_T_MASK = SF_T - 1;
  static const int SEC_PER_WEEK(604800);
  static const double EPSILON(1e-15);
  float dt;

  //Calculate IMU time
  int week_, secondsOfWeek_;
  float fractionOfSecond_;
  uint64_t tIndex = msg->tIndexTrunc;
  updateIMUtime(tIndex, week_, secondsOfWeek_, fractionOfSecond_);
  /*//setWithTruncatedSampleTime()
  int32_t tFracIndex = tIndexKconfig;
  //Modified slightly from gss->basetime.cpp
  //NOTE: truncL <=> msg->tIndexTrunc;
  //constexpr uint64_t one = 1ul; //defined at class level
  const uint64_t trunc = one << tIndex;
  const uint64_t truncHalf = trunc>>1;
  const uint64_t truncM = trunc-1;
  //ASSERT(tIndex < trunc);
  uint64_t nWholeSeconds = trefWeek*SEC_PER_WEEK + trefSecOfWeek;
  uint64_t samplesTimesDenom = nWholeSeconds*sampleFreqNum;
  float fracSecsTimesNum = trefFracSecs*sampleFreqNum;
  float fracSecsTimesNumFloor = std::floor(fracSecsTimesNum);
  samplesTimesDenom += static_cast<uint64_t>(fracSecsTimesNumFloor);
  uint64_t reftIndex = samplesTimesDenom/sampleFreqDen;
  uint64_t reftIndexRem = samplesTimesDenom%sampleFreqDen;
  float fracSamples = (reftIndexRem+(fracSecsTimesNum-fracSecsTimesNumFloor))
    /sampleFreqDen;
  int32_t reftFracIndex = static_cast<int32_t>(std::floor(fracSamples*SF_T + 0.5));
  reftIndex += (reftFracIndex >> SF_TL);
  reftFracIndex = (reftFracIndex & SF_T_MASK);
  uint64_t fulltIndex = (reftIndex & ~truncM) | tIndex;
  uint64_t truncReftIndex = reftIndex & truncM;
  if(truncReftIndex >= truncHalf){
    if(tIndex < truncReftIndex-truncHalf)
      fulltIndex += trunc;
  }
  else{
    if(tIndex > truncReftIndex+truncHalf)
      fulltIndex -= trunc;
  }

  //setWithSampleTime()
  tIndex = fulltIndex;
  const float delt = (static_cast<float>(sampleFreqDen))/sampleFreqNum;
  // One interval is equal to sampleFreqDen seconds
  const uint32_t nWholeIntervals = static_cast<uint32_t>(tIndex/sampleFreqNum);
  // Whole samples remaining after an integer number of whole intervals has
  // been removed
  const uint32_t nWholeRemainingSamples = static_cast<uint32_t>(tIndex % sampleFreqNum);
  const float secondsWithinFractionalInterval =
    (static_cast<float>(nWholeRemainingSamples) +
     (static_cast<float>(tFracIndex)/SF_T))*delt;
  const uint32_t nWholeSecondsInWholeIntervals = nWholeIntervals*sampleFreqDen;
  const uint32_t nWholeSecondsWithinFractionalInterval =
  static_cast<uint32_t>(std::floor(secondsWithinFractionalInterval));
  int week_ = nWholeSecondsInWholeIntervals/SEC_PER_WEEK;
  int secondsOfWeek_ = (nWholeSecondsInWholeIntervals % SEC_PER_WEEK) +
    nWholeSecondsWithinFractionalInterval;
  float fractionOfSecond_ = secondsWithinFractionalInterval -
    nWholeSecondsWithinFractionalInterval;

  //normalize()
  if(std::fabs(fractionOfSecond_) >= 1.0){
    int32_t sec = static_cast<int32_t>(fractionOfSecond_);
    secondsOfWeek_ += sec;
    fractionOfSecond_ -= static_cast<double>(sec);
  }
  if(std::abs(secondsOfWeek_) >= SEC_PER_WEEK){
    const uint64_t wholeWeeks = secondsOfWeek_/SEC_PER_WEEK; //const auto in normalize()
    week_ += wholeWeeks;
    secondsOfWeek_ -= (wholeWeeks*SEC_PER_WEEK);
  }
  if(std::fabs(fractionOfSecond_) < EPSILON)
    fractionOfSecond_ = 0.0;
  // Enforce non-negative fractionOfSecond_ and secondsOfWeek_; a
  // negative week_ is fine.
  while(fractionOfSecond_ < 0){
    fractionOfSecond_ += 1.0;
    --secondsOfWeek_;
  }
  while(secondsOfWeek_ < 0){
    secondsOfWeek_ += SEC_PER_WEEK;
    --week_;
  }*/

  float thisTime=week_*SEC_PER_WEEK+secondsOfWeek_+fractionOfSecond_;
  //NOTE: tLastProcessed is the last gps OR imu measurement processed whereas tLastImu is JUST imu
  std::cout << "Week: " << week_ << std::endl;
  std::cout << "Sec : " << secondsOfWeek_ << std::endl;
  std::cout << "Fsec: " << fractionOfSecond_ << std::endl;

  dt = thisTime-tLastImu;
  tLastImu=thisTime;
  if(dt<=1e-9)
  {
    //std::cout << "Error: 1ns between IMU measurements" << std::endl;
    return;
  }

  imuAccelMeas(0) = msg->acceleration[0] * imuConfigAccel;
  imuAccelMeas(1) = msg->acceleration[1] * imuConfigAccel;
  imuAccelMeas(2) = msg->acceleration[2] * imuConfigAccel;
  imuAttRateMeas(0) = msg->angularRate[0] * imuConfigAttRate;
  imuAttRateMeas(1) = msg->angularRate[1] * imuConfigAttRate;
  imuAttRateMeas(2) = msg->angularRate[2] * imuConfigAttRate;

  //Rotate gyro/accel measurements to body frame
  //NOTE1: this does not need Rpqr convention
  //NOTE2: this is hardcoded for the lynx as mounted on phoenix and company
  Eigen::Matrix3d Raccel, Rgyro;
  Raccel<<0,-1,0, -1,0,0, 0,0,-1;
  Rgyro=Raccel;
  imuAccelMeas = Raccel*imuAccelMeas;
  imuAttRateMeas = Rgyro*imuAttRateMeas;
  Eigen::Vector3d gamma0(xState(6),xState(7),xState(8));
  imuAccelMeas = imuAccelMeas - updateRBIfromGamma(RBI,gamma0)*Eigen::Vector3d(0,0,9.81);

  std::cout << "has estimated RBI: " << hasRBI <<std::endl;
  //Run CF if calibrated
  if(isCalibrated)
  {
    double dtLastProc = thisTime - tLastProcessed;
    Eigen::Matrix<double,15,15> Fmat_local = getFmatrixCF(dtLastProc,imuAccelMeas,imuAttRateMeas,RBI);
    xState=Fmat_local*xState;
    Fimu=Fmat_local*Fimu;
    RBI=updateRBIfromGamma(RBI, xState.middleRows(6,3));
    xState.middleRows(6,3)=Eigen::Vector3d::Zero();

    publishOdomAndMocap();
  }else if(hasRBI)
  {  //if RBI has been calculated but the biases have not been calculated
    ba0=ba0+1/100*(imuAccelMeas+RBI*Eigen::Vector3d(0,0,9.81)); //inefficient
    bg0=bg0+1/100*imuAttRateMeas;
    std::cout<<"calibcounter"<<counter<<std::endl;
    std::cout<<"RBI generated, estimating biases"<<std::endl;
    xState<<internal_rI(0),internal_rI(1),internal_rI(2), 0,0,0, 0,0,0, ba0(0),ba0(1),ba0(2), bg0(0),bg0(1),bg0(2);
    //std::cout<<xState<<std::endl;
    counter++;
    // Try ground calibration step for simplicity
    if(counter>=100){isCalibrated=true;}
  }

}


void gpsImuNode::PublishTransform(const geometry_msgs::Pose &pose,
                               const std_msgs::Header &header,
                               const std::string &child_frame_id)
{
  // Publish tf
  geometry_msgs::Vector3 translation;
  translation.x = pose.position.x;
  translation.y = pose.position.y;
  translation.z = pose.position.z;

  geometry_msgs::TransformStamped transform_stamped;
  transform_stamped.header = header;
  transform_stamped.child_frame_id = child_frame_id;
  transform_stamped.transform.translation = translation;
  transform_stamped.transform.rotation = pose.orientation;

  tf_broadcaster_.sendTransform(transform_stamped);
}


Eigen::Matrix3d gpsImuNode::rotMatFromEuler(Eigen::Vector3d ee)
{
  double phi=ee(0);
  double theta=ee(1);
  double psi=ee(2);
  Eigen::Matrix3d RR;
  RR<<cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta),
      sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), sin(theta)*sin(phi)*sin(psi)+cos(phi)*cos(psi), sin(phi)*cos(theta),
      cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), cos(phi)*cos(theta);
}
Eigen::Matrix3d gpsImuNode::rotMatFromQuat(Eigen::Quaterniond qq)
{
  double xx=qq.x();
  double yy=qq.y();
  double zz=qq.z();
  double ww=qq.w();
  Eigen::Matrix3d RR;
  RR << 1-2*yy*yy-2*zz*zz, 2*xx*yy+2*ww*zz, 2*xx*zz-2*ww*yy,
        2*xx*yy-2*ww*zz, 1-2*xx*xx-2*zz*zz, 2*yy*zz+2*ww*xx,
        2*xx*zz+2*ww*yy, 2*yy*zz-2*ww*xx, 1-2*xx*xx-2*yy*yy;
}


} //end namespace


int main(int argc, char **argv)
{
  ros::init(argc, argv, "gpsimu_odom");
  ros::NodeHandle nh;

  try
  {
    gpsimu_odom::gpsImuNode gpsimu_odom(nh);
    ros::spin();
  }
  catch(const std::exception &e)
  {
    ROS_ERROR("%s: %s", nh.getNamespace().c_str(), e.what());
    return 1;
  }
  return 0;
}

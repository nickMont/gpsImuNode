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
  throttleMax = tmax;

  one = 1ul;

  //Get additional parameters for the kalkman filter
  nh.param(quadName + "/max_accel", max_accel, 2.0);
  nh.param(quadName + "/publish_tf", publish_tf_, true);
  nh.param<std::string>(quadName + "/child_frame_id", child_frame_id_, "base_link");
  if(publish_tf_ && child_frame_id_.empty())
    throw std::runtime_error("gpsimu_odom: child_frame_id required for publishing tf");

  // There should only be one gps_fps, so we read from nh
  double gps_fps;
  nh.param(quadName + "/gps_fps", gps_fps, 20.0);
  ROS_ASSERT(gps_fps > 0.0);

  //should be a const but catkin doesn't like scoping it
  pi = std::atan(1.0)*4;

         /*baseECEF_vector(0) = msg->rx+msg->rxRov; //NOTE: THIS SHOULD BE READ IN VIA .LAUNCH WHEN USING GLOBAL FRAME
        baseECEF_vector(1) = msg->ry+msg->ryRov;
        baseECEF_vector(2) = msg->rz+msg->rzRov;*/
  Recef2enu=ecef2enu_rotMatrix(baseECEF_vector);
  //baseENU_vector=Recef2enu*baseECEF_vector;

  lastRTKtime=0;
  lastA2Dtime=0;
  //internalQuat.resize(4);
  internalSeq=0;
  sec_in_week = 604800;
  kfInit=false; //KF will need to be initialized
  throttleSetpoint = 9.81/throttleMax; //the floor is the throttle
  quaternionSetpoint.x()=0; quaternionSetpoint.y()=0; quaternionSetpoint.z()=0; quaternionSetpoint.w()=1;

  //verbose parameters
  ROS_INFO("max_accel: %f", max_accel);
  ROS_INFO("publish_tf_: %d", publish_tf_);
  ROS_INFO("child_frame_id: %s", child_frame_id_.c_str());
  ROS_INFO("ROS topic: %s", quadPoseTopic.c_str());
  ROS_INFO("Node name: %s", quadName.c_str());
  ROS_INFO("gps_fps: %f", gps_fps);

  // Initialize publishers and subscribers
  odom_pub_ = nh.advertise<nav_msgs::Odometry>("odom", 10); //MUST have a node namespace, ns="quadName", in launchfile
  localOdom_pub_ = nh.advertise<nav_msgs::Odometry>("local_odom", 10);
  mocap_pub_ = nh.advertise<geometry_msgs::PoseStamped>("mavros/mocap/pose", 10);
/*  gps_sub_ = nh.subscribe(quadPoseTopic, 10, &gpsImuNode::gpsCallback,
                            this, ros::TransportHints().tcpNoDelay());*/
//  internalPosePub_ = nh.advertise<geometry_msgs::PoseStamped>(posePubTopic,10);
  rtkSub_ = nh.subscribe("SingleBaselineRTK",10,&gpsImuNode::singleBaselineRTKCallback,
                            this, ros::TransportHints().tcpNoDelay());
 /* a2dSub_ = nh.subscribe("Attitude2D",10,&gpsImuNode::attitude2DCallback,
                            this, ros::TransportHints().tcpNoDelay());*/
  imuSub_ = nh.subscribe("IMU",10, &gpsImuNode::imuDataCallback,
                            this, ros::TransportHints().tcpNoDelay());
  imuConfigSub_ = nh.subscribe("IMUConfig",10, &gpsImuNode::imuConfigCallback,
                            this, ros::TransportHints().tcpNoDelay());
  tOffsetSub_ = nh.subscribe("ObservablesMeasurementTime",10,&gpsImuNode::tOffCallback,
                            this, ros::TransportHints().tcpNoDelay());
  ROS_INFO("Waiting for IMU config data, this may take a moment...");
  gbx_ros_bridge_msgs::ImuConfig::ConstPtr imuConfigMsg = ros::topic::waitForMessage<gbx_ros_bridge_msgs::ImuConfig>("IMUConfig");
  imuConfigAccel = imuConfigMsg->lsbToMetersPerSecSq;
  imuConfigAttRate = imuConfigMsg->lsbToRadPerSec;
  ROS_INFO("IMU configuration recorded, finishing startup.");
  ROS_INFO("Startup complete");

  //Get initial pose
  //initPose_ = ros::topic::waitForMessage<geometry_msgs::PoseStamped>(quadPoseTopic);
  //geometry_msgs::PoseStamped initPose_;
}


//Get reference timing from A2D since it's the most reliable
void gpsImuNode::attitude2DCallback(const gbx_ros_bridge_msgs::Attitude2D::ConstPtr &msg)
{ 
  /*
  if(msg->tSolution.week > 1) //GPS week will always be >1 if messages are being sent
  {
    trefWeek = msg->tSolution.week;
    trefSecOfWeek = msg->tSolution.secondsOfWeek;
    trefFracSecs = msg->tSolution.fractionOfSecond;
  }*/
}


void gpsImuNode::tOffCallback(const gbx_ros_bridge_msgs::ObservablesMeasurementTime::ConstPtr &msg)
{
  //tMeasOffset=msg->tOffset.week*SEC_IN_WEEK+msg->tOffset.secondsOfWeek+msg->tOffset.fractionOfSecond;
  //time offset from converting RRT to ORT
}


void gpsImuNode::imuConfigCallback(const gbx_ros_bridge_msgs::ImuConfig::ConstPtr &msg)
{
  imuConfigAccel = msg->lsbToMetersPerSecSq; //scaling to m/s2 from "non-engineering units"
  imuConfigAttRate = msg->lsbToRadPerSec; //scaling to rad/s from "non-engineering units"
  //imuSampleFreq = msg->sampleFreqNumerator/msg->sampleFreqDenominator/36/3600;  //samples per second
  
  sampleFreqNum = msg->sampleFreqNumerator;
  sampleFreqDen = msg->sampleFreqDenominator;
  
}


void gpsImuNode::imuDataCallback(const gbx_ros_bridge_msgs::Imu::ConstPtr &msg)
{
  static bool isCalibrated=false;
  static int counter=-1;
  static Eigen::Matrix<double,100,3> imuAStore, imuGStore;
  static float tLast=0;
  static Eigen::Vector3d ba0=Eigen::Vector3d(0,0,0);
  static Eigen::Vector3d bg0=Eigen::Vector3d(0,0,0);
  static Eigen::Matrix3d RBI=Eigen::Matrix3d::Identity();
  static const int SF_TL = 24;
  static const int32_t SF_T = 0x1 << SF_TL;
  static const int32_t SF_T_MASK = SF_T - 1;
  float dt;


  //const s32 SF_T = 0x1 << SF_TL;
  int32_t tFracIndex = 0; //==0 in Matthew's code
  int SEC_PER_WEEK=604800;
  //Modified slightly from gss->basetime.cpp
  //NOTE: truncL <=> msg->tIndexTrunc;
  //constexpr uint64_t one = 1ul;
  const uint64_t trunc = one << imuTimeTrunc;
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
  uint64_t tIndex = (reftIndex & ~truncM) | tIndex;
  uint64_t truncReftIndex = reftIndex & truncM;

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
  float thisTime=week_*SEC_PER_WEEK+secondsOfWeek_+fractionOfSecond_;
  dt = thisTime-tLast;

  if(dt<=1e-8)
  {
    return;
  }

  //Only update time if you accept the data
  tLast=thisTime;

  counter++;

  imuAccelMeas(0) = msg->acceleration[0] * imuConfigAccel;
  imuAccelMeas(1) = msg->acceleration[1] * imuConfigAccel;
  imuAccelMeas(2) = msg->acceleration[2] * imuConfigAccel;
  imuAttRateMeas(0) = msg->angularRate[0] * imuConfigAttRate;
  imuAttRateMeas(1) = msg->angularRate[1] * imuConfigAttRate;
  imuAttRateMeas(2) = msg->angularRate[2] * imuConfigAttRate;

  //Rotate gyro/accel to body frame
  //TODO: gyro should use Rpqr convention.
  Eigen::Matrix3d Raccel, Rgyro;
  Raccel<<-1,0,0, 0,-1,0, 0,0,-1;
  Rgyro<<1,0,0, 0,1,0, 0,0,1;
  imuAccelMeas=Raccel*imuAccelMeas;
  imuAttRateMeas=Rgyro*imuAttRateMeas;


  imuTimeTrunc = msg->tIndexTrunc;
  
  // Try ground calibration step for simplicity
  if(counter>=100){isCalibrated=true;}

  //Run CF if calibrated
  if(isCalibrated)
  {
    Eigen::Matrix<double,15,15> Fmat_local = getFmatrixCF(dt,imuAccelMeas,imuAttRateMeas,RBI);
    xState=Fmat_local*xState;
    RBI=RBI*( Eigen::Matrix3d::Identity()+hatmat(xState.middleRows(6,3)) );
    xState.middleRows(6,3)=Eigen::Vector3d::Zero();

  }else{
/*    imuAStore(counter,0)=imuAccelMeas(0);
    imuAStore(counter,1)=imuAccelMeas(1);
    imuAStore(counter,2)=imuAccelMeas(2);
    imuGStore(counter,0)=imuAttRateMeas(0);
    imuGStore(counter,1)=imuAttRateMeas(1);
    imuGStore(counter,2)=imuAttRateMeas(2);*/
    ba0=ba0+1/100*(imuAccelMeas+RBI*Eigen::Vector3d(0,0,9.81)); //inefficient
    bg0=bg0+1/100*imuAttRateMeas;

    xState<<0,0,0, 0,0,0, 0,0,0, ba0(0),ba0(1),ba0(2), bg0(0),bg0(1),bg0(2);
  }

}

/*void gpsImuNode::gpsCallback()
{
  //do GPS stuff

  //Reset for next iter
  Fimu=Eigen::Matrix<double,15,15>::Identity();
}*/


void gpsImuNode::singleBaselineRTKCallback(const gbx_ros_bridge_msgs::SingleBaselineRTK::ConstPtr &msg)
{
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

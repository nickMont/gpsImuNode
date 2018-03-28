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

  //First get rbi(0), then get biases(0)
  hasRBI = false;
  isCalibrated=false;

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
  navSub_ = nh.subscribe("NavigationSolution",10,&gpsImuNode::navsolCallback,
                            this, ros::TransportHints().tcpNoDelay());
  tOffsetSub_ = nh.subscribe("ObservablesMeasurementTime",10,&gpsImuNode::tOffsetCallback,
                            this, ros::TransportHints().tcpNoDelay());

  //Load IMU config data
  ROS_INFO("Waiting for IMU config data, this may take a moment...");
  gbx_ros_bridge_msgs::ImuConfig::ConstPtr imuConfigMsg = 
        ros::topic::waitForMessage<gbx_ros_bridge_msgs::ImuConfig>("IMUConfig");
  imuConfigAccel = imuConfigMsg->lsbToMetersPerSecSq;
  imuConfigAttRate = imuConfigMsg->lsbToRadPerSec;
  sampleFreqNum = imuConfigMsg->sampleFreqNumerator;
  sampleFreqDen = imuConfigMsg->sampleFreqDenominator;  
  tIndexConfig = imuConfigMsg->tIndexk;
  ROS_INFO("IMU configuration recorded.");

  //Load offset time data
  ROS_INFO("Waiting for offset time, this may take a moment...");
  gbx_ros_bridge_msgs::ObservablesMeasurementTime::ConstPtr toffsetMsg =
        ros::topic::waitForMessage<gbx_ros_bridge_msgs::ObservablesMeasurementTime>("ObservablesMeasurementTime");
  toffsetWeek = toffsetMsg->tOffset.week;
  toffsetSecOfWeek = toffsetMsg->tOffset.secondsOfWeek;
  toffsetFracSecs = toffsetMsg->tOffset.fractionOfSecond;
  ROS_INFO("Time offset from RRT to ORT recorded.");

  //Get dtRX0
  gbx_ros_bridge_msgs::NavigationSolution::ConstPtr navsolMsg = 
        ros::topic::waitForMessage<gbx_ros_bridge_msgs::NavigationSolution>("NavigationSolution");
  dtRX_meters = navsolMsg->deltatRxMeters;
  ROS_INFO("Time offset from RX to GPS obtained.");

  ROS_INFO("Startup complete.");
}


void gpsImuNode::navsolCallback(const gbx_ros_bridge_msgs::NavigationSolution::ConstPtr &msg)
{
  //static int SEC_PER_WEEK = 604800;

  dtRX_meters = msg->deltatRxMeters;
}


//Get reference RRT time and measurement offset time from Observables message
void gpsImuNode::tOffsetCallback(const gbx_ros_bridge_msgs::ObservablesMeasurementTime::ConstPtr &msg)
{
//  static const int SEC_PER_WEEK(604800);
//  tmeasWeek = msg->tMeasurement.week;
//  tmeasFracSecs = msg->tMeasurement.fractionOfSecond;
//  tmeasSecOfWeek = msg->tMeasurement.secondsOfWeek;

  toffsetWeek = msg->tOffset.week;
  toffsetSecOfWeek = msg->tOffset.secondsOfWeek;
  toffsetFracSecs = msg->tOffset.fractionOfSecond;

//  const double tmeas = tmeasWeek*SEC_PER_WEEK + tmeasSecOfWeek + tmeasFracSecs;
//  const double tnavsol = tnavsolWeek*SEC_PER_WEEK + tnavsolSecOfWeek + tnavsolFracSecs;
//  const double tOffset = toffsetWeek*SEC_PER_WEEK + toffsetSecOfWeek + toffsetFracSecs;
  //std::cout << "tsol_MEAS: " << tmeas+tOffset - tnavsol << std::endl;
  //std::cout << "navsol: " << tnavsolWeek << " " << tnavsolSecOfWeek << " " << tnavsolFracSecs <<std::endl;
  //std::cout << "calcd : " << tmeasWeek+toffsetWeek << " " << tmeasSecOfWeek+toffsetSecOfWeek << " " << tmeasFracSecs+toffsetFracSecs <<std::endl;
}


//Get upper 32 bits of tIndex counter
void gpsImuNode::imuConfigCallback(const gbx_ros_bridge_msgs::ImuConfig::ConstPtr &msg)
{
  ROS_INFO("Config message received.");
  imuConfigAccel = msg->lsbToMetersPerSecSq; //scaling to m/s2 from "non-engineering units"
  imuConfigAttRate = msg->lsbToRadPerSec; //scaling to rad/s from "non-engineering units"
  //imuSampleFreq = msg->sampleFreqNumerator/msg->sampleFreqDenominator/36/3600;  //samples per second
  
  sampleFreqNum = msg->sampleFreqNumerator;
  sampleFreqDen = msg->sampleFreqDenominator;
  tIndexConfig = msg->tIndexk;
}


void gpsImuNode::imuDataCallback(const gbx_ros_bridge_msgs::Imu::ConstPtr &msg)
{
  //Initialization variables
  static int counter=0;
  static Eigen::Matrix<double,100,3> imuAStore=Eigen::MatrixXd::Zero(100,3);
  static Eigen::Matrix<double,100,3> imuGStore=Eigen::MatrixXd::Zero(100,3);
  static float tLastImu=0;
  static Eigen::Vector3d ba0=Eigen::Vector3d(0,0,0);
  static Eigen::Vector3d bg0=Eigen::Vector3d(0,0,0);
  //gps variables
  static const int SEC_PER_WEEK(604800);
  static const double cLight(299792458);
  float dt;

  //Calculate IMU time
  //int week_, secondsOfWeek_;
  //float fractionOfSecond_;
  uint64_t tIndex = msg->tIndexTrunc;
  //updateIMUtimeRRT(tIndex, week_, secondsOfWeek_, fractionOfSecond_);
  //const float thisTimeRRT = week_*SEC_PER_WEEK+secondsOfWeek_+fractionOfSecond_;
  //uint64_t tFullIndex = tIndexConfig << 32 + tIndex;
  //float thisTimeRRT = tFullIndex*sampleFreqNum/sampleFreqDen;
  //std::cout << thisTimeRRT;
  //const float thisTimeORT = thisTimeRRT + toffsetWeek*SEC_PER_WEEK + toffsetSecOfWeek + toffsetFracSecs;
  //const float thisTime = thisTimeORT - deltaRX/cLight;
  long long int mask = 0xffffffff; // This is: (1 << 32) - 1
  long long int tIndexFull = (tIndexConfig & ~mask) | (tIndex & mask); // You don't want to bit-shift by 32-bits. You want to bit-mask to use the upper 32-bits from tIndexImuConfig and the lower 32-bits from tIndexImu
  double sampleFreq = sampleFreqNum/sampleFreqDen;
  double tRRT = tIndexFull/sampleFreq; // tIndexFull is in units of samples. Divide by sample rate to get units of seconds.
  //std::cout << "rrt seconds:" << tRRT << std::endl;
  //std::cout << "meas-rrt:   " << tRRT - (tmeasSecOfWeek + tmeasFracSecs) << std::endl;
  //const double tmeas = tmeasFracSecs + tmeasWeek*SEC_PER_WEEK + tmeasSecOfWeek;
  //const double tnavsol = tnavsolFracSecs + tnavsolWeek*SEC_PER_WEEK + tnavsolSecOfWeek;
  const double tOffset = toffsetFracSecs + toffsetWeek*SEC_PER_WEEK + toffsetSecOfWeek;
  const double tORT = tRRT + tOffset; // tOffset comes from ObservablesMeasurementTime
  const double tGPS = tORT - dtRX_meters/cLight; // dtrx comes from NavigationSolution
  const double thisTime = tGPS;
 
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
    //std::cout<<"RBI generated, estimating biases"<<std::endl;
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

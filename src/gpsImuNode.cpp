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
  Qimu=1.0e-4*Eigen::Matrix<double,6,6>::Identity();
  Rk=1.0e-3*Eigen::Matrix<double,6,6>::Identity();

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
  Eigen::Matrix<double,15,1> PdiagElements;
  PdiagElements << 1.0e-2,1.0e-2,1.0e-2, 1.0e-3,1.0e-3,1.0e-3,
      1.0e-2,1.0e-2,1.0e-2, 1.0e-1,1.0e-1,1.0e-1, 1.0e-1, 1.0e-1, 1.0e-1;
  Pimu = PdiagElements.asDiagonal();
  R_G2wrw=Rwrw*Recef2enu;
  RBI=Eigen::MatrixXd::Identity(3,3);

  //THIS SHOULD BE INITIALIZED IN GPSCALLBACKS
  rRefImu<<0,0,0; //location of rI_0 for imu

  //Stuff for rI
  sec_in_week = 604800;
  lastRTKtime=0;
  lastA2Dtime=0;
  internalSeq=0;

  //Secondary to primary vector in body frame
  l_s2p<<0.195,0,0;
  l_imu<<-0.0975,0,-0.05;
  l_cg2p<<0.0975,0,0.05;

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


void gpsImuNode::imuDataCallback(const gbx_ros_bridge_msgs::Imu::ConstPtr &msg)
{
  //Initialization variables
  static int counter=0;
  static uint64_t imuSeq=0;
  static Eigen::Matrix<double,100,3> imuAStore=Eigen::MatrixXd::Zero(100,3);
  static Eigen::Matrix<double,100,3> imuGStore=Eigen::MatrixXd::Zero(100,3);
  static double tLastImu=0;
  static Eigen::Vector3d ba0=Eigen::Vector3d(0,0,0);
  static Eigen::Vector3d bg0=Eigen::Vector3d(0,0,0);
  //gps variables
  static const int SEC_PER_WEEK(604800);
  static const double cLight(299792458);
  static const long long int mask = 0xffffffff; // This is: (1 << 32) - 1
  double dt;

  //Calculate IMU time
  const uint64_t tIndex = msg->tIndexTrunc;
  const long long int tIndexFull = (tIndexConfig & ~mask) | (tIndex & mask); //Bit-mask, don't bit-shift
  const double sampleFreq = sampleFreqNum/sampleFreqDen;
  const double tRRT = tIndexFull/sampleFreq; //tIndexFull is in samples, divide by samples/s
  const double tOffset = toffsetFracSecs + toffsetWeek*SEC_PER_WEEK + toffsetSecOfWeek;
  const double tORT = tRRT + tOffset; //tOffset comes from ObservablesMeasurementTime
  const double tGPS = tORT - dtRX_meters/cLight; //dtrx comes from NavigationSolution
  const double thisTime = tGPS;
  std::cout.precision(17);
  dt = thisTime-tLastImu;
  if(dt<=1e-9)
  {
    std::cout << "Error: 1ns between IMU measurements" << std::endl;
    //std::cout << "Time: " << thisTime-double(toffsetWeek*SEC_PER_WEEK) <<std::endl;
    return;
  }
  //Only update last time used IF this time is accepted
  tLastImu=thisTime;

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
  //attRateMeasOrig = imuAttRateMeas;
  //accelMeasOrig = imuAccelMeas;
  //Eigen::Vector3d correctedImuAccelMeas = imuAccelMeas - updateRBIfromGamma(RBI,xState.middleRows(6,3))*Eigen::Vector3d(0,0,9.81);

  //Run CF if calibrated
  if(isCalibrated)
  {
    imuSeq++;
    double dtLastProc = thisTime - tLastProcessed;
    if(dtLastProc>0) //Just in case a gps message is received late
    { 
      //propagate state nonlinearly
      xState=fdyn(xState,dtLastProc,imuAccelMeas,imuAttRateMeas,RBI,l_imu);
      RBI=updateRBIfromGamma(RBI, xState.middleRows(6,3));
      xState.middleRows(6,3)=Eigen::Vector3d::Zero();
      //Augment F matrix
      Eigen::Matrix<double,15,15> Fmat_local = getFmatrixCF(dtLastProc,imuAccelMeas,imuAttRateMeas,RBI);
      //Fmat_local = getNumderivF(double(1e-9), dtLastProc, xState, accelMeasOrig, attRateMeasOrig,RBI, l_imu);
      Fimu=Fmat_local*Fimu;
      std::cout<<"rbidet: " << RBI.determinant() << std::endl;
      std::cout << "RBI:" << std::endl << RBI <<std::endl;
      //std::cout<<"xk"<<std::endl<<xState.topRows(3)<<std::endl;
      publishOdomAndMocap();
      tLastProcessed = thisTime;
      counter++;
      if(counter%800==0)
      {
        RBI = orthonormalize(RBI);
      }
    }
  }else if(hasRBI)
  { 
    //std::cout << "bg" <<std::endl<<bg0 <<std::endl;
    //if RBI has been calculated but the biases have not been calculated
    ba0=ba0+0.01*(imuAccelMeas - RBI*Eigen::Vector3d(0,0,9.81)); //inefficient
    bg0=bg0+0.01*imuAttRateMeas;
    //std::cout<<"RBI generated, estimating biases"<<std::endl;
    Eigen::Vector3d rI0;
    rI0 = internal_rI - RBI.transpose()*l_cg2p;
    xState<<rI0(0),rI0(1),rI0(2), 0,0,0, 0,0,0, ba0(0),ba0(1),ba0(2), bg0(0),bg0(1),bg0(2);
    //std::cout<<xState<<std::endl;
    counter++;
    // Try ground calibration step for simplicity
    if(counter>=100)
    {
      isCalibrated = true;
      tLastProcessed = thisTime;
      std::cout << "state0" <<std::endl<<xState<<std::endl;
    }
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

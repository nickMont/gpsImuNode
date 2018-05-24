#include <Eigen/Geometry>
#include "estimationNode.hpp"
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
estimationNode::estimationNode(ros::NodeHandle &nh)
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

  //Covariances
  tauA=1000.0;
  tauG=1000.0;
  double pA,pG;
  pA=1.0e-2;
  pG=1.0e-3;
  Eigen::Matrix<double,6,6> temp;
  Eigen::Matrix<double,6,1> P6diagElements;
  P6diagElements << pA*0.000241,pA*0.000241,pA*0.000241, pG*3.0e-3,pG*3.0e-3,pG*3.0e-3;
  Qimu = P6diagElements.asDiagonal();
  P6diagElements << 0.000036,0.000036,0.000144, 0.000144,0.000144,0.000144;  //Params file
  P6diagElements << 0.000036,0.000036,0.000019, 0.000144,0.000144,0.000144;  //Test values
  Rk = P6diagElements.asDiagonal();
  Qk12 = Eigen::Matrix<double,12,12>::Zero();

  /*//gpsMeas class tester
  Eigen::Vector3d testvec(1.0,2.0,3.0);
  gpsMeas g2;
  g2.setMeas(0.0,&testvec,&testvec);
  Eigen::Vector3d out1,out2;
  double t2;
  g2.getMeas(t2,out1,out2);
  std::cout<<t2<<std::endl<<out1<<std::endl<<out2<<std::endl;*/
    

  //12x12 Q covariance
  const double dtIMU=1.0/73.25; //from (rostopic hz /phoenix/imu -r 10)/10
  const double alphaA = exp(-dtIMU/tauA);
  const double gScale = 9.81/1000.0;
  const double thetaScale = pi/180.0;
  const double alphaG = exp(-dtIMU/tauG);
  //Covariance elements: gyro output, gyro bias, accel output, accel bias

  //Filled spec sheet data
  Qk12.topLeftCorner(3,3) = pow(0.1*pi/180/sqrt(dtIMU),2)*Eigen::Matrix3d::Identity(); //see datasheet  
  Qk12.block(3,3,3,3) = pow(thetaScale*100.0/360.0,2)*(1.-alphaG*alphaG)*Eigen::Matrix3d::Identity(); //random tests
  Qk12.block(6,6,3,3) = pow(9.81/1.0e6*150.0/sqrt(dtIMU),2)*Eigen::Matrix3d::Identity(); //see datasheet
  Qk12.bottomRightCorner(3,3) = pow(gScale*1000.0,2)*(1.-alphaA*alphaA)*Eigen::Matrix3d::Identity(); //random tests

  //Testing based on ground test.  Works but velocity is noisy
  Qk12.topLeftCorner(3,3) = 6.95e-4*Eigen::Matrix3d::Identity();
  Qk12.block(3,3,3,3) = pow(thetaScale*100.0/360.0,2)*(1.-alphaG*alphaG)*Eigen::Matrix3d::Identity();
  Qk12.block(6,6,3,3) = 0.0045*Eigen::Matrix3d::Identity();
  Qk12.bottomRightCorner(3,3) = 0.1*pow(gScale*1000.0,2)*(1.-alphaA*alphaA)*Eigen::Matrix3d::Identity();

  //Testing
  Eigen::Matrix3d QangularAccel = 0.1*Eigen::Matrix3d::Identity(); //Q caused by angular acceleration
  Eigen::Matrix3d QgyroOutput2 = 1e-3*Eigen::Matrix3d::Identity();
  Qk12.topLeftCorner(3,3) = 6.95e-4*Eigen::Matrix3d::Identity() + QgyroOutput2;
  Qk12.block(3,3,3,3) = 1.0e-6*pow(thetaScale*100.0/360.0,2)*(1.-alphaG*alphaG)*Eigen::Matrix3d::Identity();
  Qk12.block(6,6,3,3) = 0.0045*Eigen::Matrix3d::Identity() + QangularAccel;
  Qk12.bottomRightCorner(3,3) = 1.0e-6*pow(gScale*80.0,2)*(1.-alphaA*alphaA)*Eigen::Matrix3d::Identity();

  //Testing, not currently in use
  Qk12dividedByDt = Qk12/dtIMU;  //When used in filter, multiply by dt of each update step

  Pimu=Eigen::MatrixXd::Identity(15,15);
  Eigen::Matrix<double,15,1> P15diagElements;
  P15diagElements << 1.0e-4,1.0e-4,1.0e-4, 1.0e-6,1.0e-6,1.0e-6,
      1.0e-5,1.0e-5,1.0e-5, 1.0e-5,1.0e-5,1.0e-5, 1.0e-8, 1.0e-8, 1.0e-8;
  Pimu = P15diagElements.asDiagonal();
  P_report=Pimu;
  R_G2wrw=Rwrw*Recef2enu;
  RBI=Eigen::MatrixXd::Identity(3,3);

  //THIS SHOULD BE INITIALIZED IN GPSCALLBACKS
  rRefImu<<0,0,0; //location of rI_0 for imu

  //Stuff for rI
  sec_in_week = 604800.0;
  lastRTKtime=0;
  lastA2Dtime=0;
  internalSeq=0;

  //Secondary to primary vector in body frame
  l_s2p<< 0.190,0,0;
  l_imu<< -0.0884,0.0134,-0.0399;
  l_cg2p<< 0.1013,-0.0004,0.0472;

  //First get rbi(0), then get biases(0)
  hasRBI = false;
  isCalibrated=false;

  //Try to use buffer to fix timing issues. Experimental
  tryToUseBuffer = false;

  // Initialize publishers and subscribers
  //odom_pub_ = nh.advertise<nav_msgs::Odometry>("odom", 10); //MUST have a node namespace, ns="quadName", in launchfile
  localOdom_pub_ = nh.advertise<nav_msgs::Odometry>("local_odom_INS", 10);
  mocap_pub_ = nh.advertise<geometry_msgs::PoseStamped>("mavros/mocap/pose_INS", 10);
/*  gps_sub_ = nh.subscribe(quadPoseTopic, 10, &estimationNode::gpsCallback,
                            this, ros::TransportHints().tcpNoDelay());*/
//  internalPosePub_ = nh.advertise<geometry_msgs::PoseStamped>(posePubTopic,10);
  rtkSub_ = nh.subscribe("SingleBaselineRTK",10,&estimationNode::singleBaselineRTKCallback,
                            this, ros::TransportHints().tcpNoDelay());
  a2dSub_ = nh.subscribe("Attitude2D",10,&estimationNode::attitude2DCallback,
                            this, ros::TransportHints().tcpNoDelay());
  imuSub_ = nh.subscribe("IMU",10, &estimationNode::imuDataCallback,
                            this, ros::TransportHints().tcpNoDelay());
  imuConfigSub_ = nh.subscribe("IMUConfig",10, &estimationNode::imuConfigCallback,
                            this, ros::TransportHints().tcpNoDelay());
  navSub_ = nh.subscribe("NavigationSolution",10,&estimationNode::navsolCallback,
                            this, ros::TransportHints().tcpNoDelay());
  tOffsetSub_ = nh.subscribe("ObservablesMeasurementTime",10,&estimationNode::tOffsetCallback,
                            this, ros::TransportHints().tcpNoDelay());

  //Load IMU config data, establish saturations
  ROS_INFO("Waiting for IMU config data, this may take a moment...");
  gbx_ros_bridge_msgs::ImuConfig::ConstPtr imuConfigMsg = 
        ros::topic::waitForMessage<gbx_ros_bridge_msgs::ImuConfig>("IMUConfig");
  imuConfigAccel = imuConfigMsg->lsbToMetersPerSecSq;
  imuConfigAttRate = imuConfigMsg->lsbToRadPerSec;
  sampleFreqNum = imuConfigMsg->sampleFreqNumerator;
  sampleFreqDen = imuConfigMsg->sampleFreqDenominator;  
  tIndexConfig = imuConfigMsg->tIndexk;
  //maxBa = imuConfigAccel * 250.0; //imu report units times scalefactor
  //maxBg = imuConfigAttRate * 250.0;
  maxBa = imuConfigAccel*1.0e4;
  maxBg = imuConfigAttRate*1.0e4;
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


void estimationNode::imuDataCallback(const gbx_ros_bridge_msgs::Imu::ConstPtr &msg)
{
  //Initialization variables
  static int counter=0;
  static uint64_t imuSeq=0;
  static double tLastImu=0;
  static Eigen::Vector3d ba0=Eigen::Vector3d(0,0,0);
  static Eigen::Vector3d bg0=Eigen::Vector3d(0,0,0);
  //gps variables
  static const int SEC_PER_WEEK(604800);
  static const double cLight(299792458);
  static const long long int mask = 0xffffffff; // This is: (1 << 32) - 1
  static int warnCounter=0;
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
  if(dt<=1e-9) //Note: NOT abs(dt), as this checks whether or not messages are received out of order as well
  {
    std::cout << "Error: 1ns between IMU measurements" << std::endl;
    //std::cout << "Time: " << thisTime-double(toffsetWeek*SEC_PER_WEEK) <<std::endl;
    return;
  }
  //Only update last time used IF this time is accepted
  tLastImu=thisTime;

  //Use scale parameter
  imuAccelMeas(0) = msg->acceleration[0] * imuConfigAccel;
  imuAccelMeas(1) = msg->acceleration[1] * imuConfigAccel;
  imuAccelMeas(2) = msg->acceleration[2] * imuConfigAccel;
  imuAttRateMeas(0) = msg->angularRate[0] * imuConfigAttRate;
  imuAttRateMeas(1) = msg->angularRate[1] * imuConfigAttRate;
  imuAttRateMeas(2) = msg->angularRate[2] * imuConfigAttRate;

  //Rotate gyro/accel measurements to body frame
  //NOTE1: this does not need Rpqr convention
  //NOTE2: this is hardcoded for the lynx as mounted on phoenix and company
  imuAccelPrev = imuAccelMeas;
  imuAttRatePrev = imuAttRateMeas;
  Eigen::Matrix3d Raccel, Rgyro;
  Raccel << 0,-1,0, -1,0,0, 0,0,-1;
  Rgyro  << 0,-1,0, -1,0,0, 0,0,-1;
  imuAccelMeas = Raccel*imuAccelMeas;
  imuAttRateMeas = Rgyro*imuAttRateMeas;

  //Run CF if calibrated
  if(isCalibrated)
  {
    imuSeq++;
    double dtLastProc = thisTime - tLastProcessed;
    if(dtLastProc>0) //Just in case a gps message is received late
    { 
      //Propagate state nonlinearly
      /*xState=fdyn(xState,dtLastProc,imuAccelMeas,imuAttRateMeas,RBI,l_imu);
      
      //Augment F matrix
      //Eigen::Matrix<double,15,15> Fmat_local = getFmatrixCF(dtLastProc,imuAccelMeas,imuAttRateMeas,RBI);
      Eigen::Matrix<double,15,15> Fmat_local = getFmatrixCF(dtLastProc,imuAccelMeas,imuAttRateMeas,RBI,xState,l_imu);
      //Fmat_local = getNumderivF(double(1e-9), dtLastProc, xState, accelMeasOrig, attRateMeasOrig,RBI, l_imu);
      Fimu=Fmat_local*Fimu;
      
      //Calculate reported covariance
      Eigen::Matrix<double,15,6> gammak=getGammakmatrixCF(dtLastProc,RBI); //dt for gammak is from last gps to current gps
      P_report = Fmat_local*P_report*Fmat_local.transpose() + gammak*Qimu*gammak.transpose();*/

      //if(imuAccelMeas.norm()>9.81*2.0)
      //  {std::cout << "Multiple gs acceleration: " << imuAccelMeas.norm()/9.81 << std::endl;}

      //Data storage for poor man's buffer
      xStatePrev = xState;
      PimuPrev = Pimu;
      tProcPrev = tLastProcessed;

      //Propagate
      spkfPropagate15(xState,Pimu,Qk12,dtLastProc,imuAccelMeas,imuAttRateMeas,RBI,l_imu,Pimu,xState);
      //Publish
      updateType = "imu";
      P_report = Pimu;
      publishOdomAndMocap();

      //Cleanup
      tLastProcessed = thisTime;

      //Orthonormalize every ~10s
      counter++;
      /*if(counter%2==0)
      {
        RBI = orthonormalize(RBI);
      }*/

      //Warn if signal is lost
      if( (thisTime-lastRTKtime > 0.5) || (thisTime-lastA2Dtime > 0.5) )
      {
        warnCounter++;
        if(warnCounter%40==0)
        {ROS_INFO("GPS outage warning!");}
      }else
      {
        if(warnCounter!=0)
        {ROS_INFO("GPS restored.");}
        warnCounter=0;
      }
    }
  }else if(hasRBI)  //if RBI has been calculated but the biases have not been calculated
  { 
    ba0=ba0+0.01*(imuAccelMeas - RBI*Eigen::Vector3d(0,0,9.8)); //inefficient
    bg0=bg0+0.01*imuAttRateMeas;
    Eigen::Vector3d rI0;
    rI0 = internal_rI - RBI.transpose()*l_cg2p;
    initBA = ba0;
    initBG = bg0;
    xState<<rI0(0),rI0(1),rI0(2), 0,0,0, 0,0,0, ba0(0),ba0(1),ba0(2), bg0(0),bg0(1),bg0(2);
    counter++;
    // Try ground calibration step for simplicity
    if(counter>=100)
    {
      isCalibrated = true;
      tLastProcessed = thisTime;
      //std::cout << "state0" <<std::endl<<xState<<std::endl;
    }
  }

}


void estimationNode::PublishTransform(const geometry_msgs::Pose &pose,
                               const std_msgs::Header &header,
                               const std::string &child_frame_id)
{
  //Publish tf
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


} //end namespace


int main(int argc, char **argv)
{
  ros::init(argc, argv, "gpsimu_odom");
  ros::NodeHandle nh;

  try
  {
    gpsimu_odom::estimationNode gpsimu_odom(nh);
    ros::spin();
  }
  catch(const std::exception &e)
  {
    ROS_ERROR("%s: %s", nh.getNamespace().c_str(), e.what());
    return 1;
  }
  return 0;
}

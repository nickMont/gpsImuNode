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
    ros::param::get(quadName + "/arenaCenterX", zeroInECEF_(0));
    ros::param::get(quadName + "/arenaCenterY", zeroInECEF_(1));
    ros::param::get(quadName + "/arenaCenterZ", zeroInECEF_(2));
    ros::param::get(quadName + "/arenaCenterX_ENU", zeroInWRW_(0));
    ros::param::get(quadName + "/arenaCenterY_ENU", zeroInWRW_(1));
    ros::param::get(quadName + "/arenaCenterZ_ENU", zeroInWRW_(2));
    ros::param::get(quadName + "/rtktopic", rtktopic);
    ros::param::get(quadName + "/a2dtopic", a2dtopic);
    ros::param::get(quadName + "/posePubTopic", posePubTopic);
    ros::param::get(quadName + "/minimumTestStat",minTestStat_);
    ros::param::get(quadName + "/maxThrust",tmax);

    //Get additional parameters for the kalkman filter
    nh.param(quadName + "/publish_tf", publish_tf_, true);
    nh.param<std::string>(quadName + "/child_frame_id", child_frame_id_, "base_link");
    if(publish_tf_ && child_frame_id_.empty())
        throw std::runtime_error("gpsimu_odom: child_frame_id required for publishing tf");

    //Should be a const but catkin doesn't like scoping it
    pi = std::atan(1.0)*4;

    Recef2enu_=ecef2enu_rotMatrix(zeroInECEF_);

    //Account for WRW rotation wrt ENU
    double thetaWRW;
    thetaWRW = 6.2*pi/180; //angle of rooftop coordinate system WRT ENU
    Rwrw_ << cos(thetaWRW), -1*sin(thetaWRW), 0,
                    sin(thetaWRW), cos(thetaWRW), 0,
                    0, 0, 1;
    Recef2wrw_ = Rwrw_*Recef2enu_;
    Recef2wrw_=Rwrw_*Recef2enu_;
    RBI_=Eigen::MatrixXd::Identity(3,3);

    //Delcarations before filling in covs
    Eigen::Matrix<double,6,6> Rk;
    Eigen::Matrix<double,12,12> Qk12, Qk12dividedByDt;
    Eigen::Matrix<double,15,15> Pimu;
    Eigen::Matrix<double,6,1> P6diagElements;
    Eigen::Matrix<double,15,1> P15diagElements;
    double tauA, tauG, maxBa, maxBg;
    Pimu=Eigen::Matrix<double,15,15>::Identity();

    //Covariances
    tauA=1000.0;
    tauG=1000.0;    
    //P6diagElements << 0.000036,0.000036,0.000144, 0.000144,0.000144,0.000144;  //Params file
    P6diagElements << 0.000036,0.000036,0.000019, 0.000144,0.000144,0.000144;  //Test values
    Rk = P6diagElements.asDiagonal();
        

    //IMU params
    const double dtIMU=1.0/73.25; //from (rostopic hz /phoenix/imu -r 10)/10
    const double alphaA = exp(-dtIMU/tauA);
    const double gScale = 9.81/1000.0;
    const double thetaScale = pi/180.0;
    const double alphaG = exp(-dtIMU/tauG);
    //Covariance elements: gyro output, gyro bias, accel output, accel bias

/*  //Filled spec sheet data
    Qk12.topLeftCorner(3,3) = pow(0.1*pi/180/sqrt(dtIMU),2)*Eigen::Matrix3d::Identity(); //see datasheet  
    Qk12.block(3,3,3,3) = pow(thetaScale*100.0/360.0,2)*(1.-alphaG*alphaG)*Eigen::Matrix3d::Identity(); //random tests
    Qk12.block(6,6,3,3) = pow(9.81/1.0e6*150.0/sqrt(dtIMU),2)*Eigen::Matrix3d::Identity(); //see datasheet
    Qk12.bottomRightCorner(3,3) = pow(gScale*1000.0,2)*(1.-alphaA*alphaA)*Eigen::Matrix3d::Identity(); //random tests */

/*  //Testing based on ground test.  Works but velocity is noisy
    Qk12.topLeftCorner(3,3) = 6.95e-4*Eigen::Matrix3d::Identity();
    Qk12.block(3,3,3,3) = pow(thetaScale*100.0/360.0,2)*(1.-alphaG*alphaG)*Eigen::Matrix3d::Identity();
    Qk12.block(6,6,3,3) = 0.0045*Eigen::Matrix3d::Identity();
    Qk12.bottomRightCorner(3,3) = 0.1*pow(gScale*1000.0,2)*(1.-alphaA*alphaA)*Eigen::Matrix3d::Identity(); */

    //Best tracking from ground set
    Eigen::Matrix3d QangularAccel = 0.1*Eigen::Matrix3d::Identity(); //Q caused by angular acceleration
    Eigen::Matrix3d QgyroOutput2 = 1e-3*Eigen::Matrix3d::Identity();
    Qk12=Eigen::Matrix<double,12,12>::Zero();
    Qk12.topLeftCorner(3,3) = 6.95e-4*Eigen::Matrix3d::Identity() + QgyroOutput2;
    Qk12.block(3,3,3,3) = 1.0e-6*pow(thetaScale*100.0/360.0,2)*(1.-alphaG*alphaG)*Eigen::Matrix3d::Identity();
    Qk12.block(6,6,3,3) = 0.0045*Eigen::Matrix3d::Identity() + QangularAccel;
    Qk12.bottomRightCorner(3,3) = 1.0e-6*pow(gScale*80.0,2)*(1.-alphaA*alphaA)*Eigen::Matrix3d::Identity();
    Qk12dividedByDt = Qk12/dtIMU;  //When used in filter, multiply by dt of each update step

    QgyroOutput_=QgyroOutput2;
    
    P15diagElements << 1.0e-4,1.0e-4,1.0e-4, 1.0e-6,1.0e-6,1.0e-6,
            1.0e-5,1.0e-5,1.0e-5, 1.0e-5,1.0e-5,1.0e-5, 1.0e-8, 1.0e-8, 1.0e-8;
    Pimu=P15diagElements.asDiagonal();

    //Stuff for rI
    sec_in_week_ = 604800.0;
    lastRTKtime_=0;
    lastA2Dtime_=0;
    internalSeq=0;

    //Secondary to primary vector in body frame
    Ls2p_<< 0.190,0,0;
    Lcg2imu_<< -0.0884,0.0134,-0.0399;
    Lcg2p_<< 0.1013,-0.0004,0.0472; //kept as class var for use in 

    //First get rbi(0), then get biases(0)
    rbiIsInitialized_ = false;
    isCalibrated_=false;

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
    imuConfigAccel_ = imuConfigMsg->lsbToMetersPerSecSq;
    imuConfigAttRate_ = imuConfigMsg->lsbToRadPerSec;
    sampleFreqNum_ = imuConfigMsg->sampleFreqNumerator;
    sampleFreqDen_ = imuConfigMsg->sampleFreqDenominator;  
    tIndexConfig_ = imuConfigMsg->tIndexk;
    //maxBa = imuConfigAccel_ * 250.0; //imu report units times scalefactor
    //maxBg = imuConfigAttRate_ * 250.0;
    maxBa = imuConfigAccel_*1.0e4;
    maxBg = imuConfigAttRate_*1.0e4;
    ROS_INFO("IMU configuration recorded.");

    //Load offset time data
    ROS_INFO("Waiting for offset time, this may take a moment...");
    gbx_ros_bridge_msgs::ObservablesMeasurementTime::ConstPtr toffsetMsg =
                ros::topic::waitForMessage<gbx_ros_bridge_msgs::ObservablesMeasurementTime>("ObservablesMeasurementTime");
    toffsetWeek_ = toffsetMsg->tOffset.week;
    toffsetSecOfWeek_ = toffsetMsg->tOffset.secondsOfWeek;
    toffsetFracSecs_ = toffsetMsg->tOffset.fractionOfSecond;
    ROS_INFO("Time offset from RRT to ORT recorded.");

    //Get dtRX0
    gbx_ros_bridge_msgs::NavigationSolution::ConstPtr navsolMsg = 
                ros::topic::waitForMessage<gbx_ros_bridge_msgs::NavigationSolution>("NavigationSolution");
    dtRXinMeters_ = navsolMsg->deltatRxMeters;
    ROS_INFO("Time offset from RX to GPS obtained.");

    ROS_INFO("Setting hard parameters for complementary filter.");
    imuFilter_.setLevers(Lcg2p_,Lcg2imu_,Ls2p_);
    imuFilter_.setCovariances(dtIMU, Qk12, Pimu, Rk);
    imuFilter_.setImuParams(tauA,tauG);
    imuFilter_.setBiasSaturationLimits(maxBa, maxBg);

    ROS_INFO("Startup complete.");
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

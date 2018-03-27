#include <Eigen/Geometry>
#include "gpsImuNode.hpp"
#include <string>
#include <iostream>
//Contains the two "wordy" A2D and SBRTK callbacks.
//Also contains some rotation matrix construction

namespace gpsimu_odom
{
void gpsImuNode::singleBaselineRTKCallback(const gbx_ros_bridge_msgs::SingleBaselineRTK::ConstPtr &msg)
{
    double ttime=msg->tSolution.secondsOfWeek + msg->tSolution.fractionOfSecond
        + msg->tSolution.week * sec_in_week -msg->deltRSec;
    if(ttime>lastRTKtime)  //only use newest time
    {
        hasAlreadyReceivedRTK=true;
        lastRTKtime=ttime;
        if(msg->testStat > minTestStat)
        {
            validRTKtest=true;
            Eigen::Vector3d tmpvec;
            //Rotate rECEF to rI and store in internal_rI
            tmpvec(0) = msg->rx + msg->rxRov - baseECEF_vector(0); //error vector from ECEF at init time
            tmpvec(1) = msg->ry + msg->ryRov - baseECEF_vector(1);
            tmpvec(2) = msg->rz + msg->rzRov - baseECEF_vector(2);
            internal_rI = Rwrw*(0.5*Recef2enu*tmpvec - n_err);
            //ROS_INFO("%f %f %f %f",msg->rx, msg->rxRov, tmpvec(0), internalPose(0)); //debugging
        }else{validRTKtest=false;}

        //If the time is new, both messages for this time have been received, and teststats are good
        if(abs(lastRTKtime-lastA2Dtime)<.001 && validA2Dtest && validRTKtest
                && hasAlreadyReceivedRTK && hasAlreadyReceivedA2D)  //only resend pose if new
        {
            internalSeq++;

            //Run CF
            Fimu = getFmatrixCF(ttime-tLastProcessed,imuAccelMeas,imuAttRateMeas,RBI) * Fimu;
            runCF(ttime-tLastProcessed);

            //Publish messages
            publishOdomAndMocap();

            //Clean up
            Fimu = Eigen::MatrixXd::Identity(15,15);
            tLastProcessed = ttime;
            RBI = updateRBIfromGamma(RBI,xState.middleRows(6,3));
            xState.middleRows(6,3)=Eigen::Vector3d::Zero();
            rRefImu=rRefImu+xState.topRows(3);
            xState.topRows(3)=Eigen::Vector3d::Zero();

            //Reset to avoid publishing twice
            hasAlreadyReceivedRTK=false; hasAlreadyReceivedA2D=false; 

        }else{ lastRTKtime=msg->tSolution.secondsOfWeek+msg->tSolution.fractionOfSecond-msg->deltRSec +
                        msg->tSolution.week * sec_in_week;}
    }
}


void gpsImuNode::attitude2DCallback(const gbx_ros_bridge_msgs::Attitude2D::ConstPtr &msg)
{
    static int rCCalibCounter=0;
    static int calibSamples=20;

      //Ignore zero messages
    if(msg->tSolution.week<=1)
    {return;}

    if(!hasRBI && msg->testStat>=100)
    {
        Eigen::Vector3d constrainedBaselineECEF(msg->rx,msg->ry,msg->rz);
        Eigen::Vector3d constrainedBaselineENU = Recef2enu*constrainedBaselineECEF;
        Eigen::Vector3d constrainedBaselineI = unit3(Rwrw*constrainedBaselineENU);
        rCtildeCalib(rCCalibCounter%calibSamples,0)=constrainedBaselineI(0);
        rCtildeCalib(rCCalibCounter%calibSamples,1)=constrainedBaselineI(1);
        rCtildeCalib(rCCalibCounter%calibSamples,2)=constrainedBaselineI(2);
        rBCalib(rCCalibCounter%calibSamples,0)=1;
        rBCalib(rCCalibCounter%calibSamples,1)=0;
        rBCalib(rCCalibCounter%calibSamples,2)=0;
        //rCB=???
        rCCalibCounter++;

        //estimate initial RBI from wahba
        if(rCCalibCounter>=calibSamples)
        {
            rCtildeCalib(calibSamples,0)=0; rBCalib(calibSamples,0)=0;
            rCtildeCalib(calibSamples,1)=0; rBCalib(calibSamples,1)=0;
            rCtildeCalib(calibSamples,2)=1; rBCalib(calibSamples,2)=1;
            Eigen::MatrixXd weights;
            weights.resize(calibSamples+1,1);
            weights.topRows(calibSamples)=0.5*1/calibSamples*Eigen::MatrixXd::Ones(calibSamples,1);
            weights(calibSamples)=0.5;
            RBI=rotMatFromWahba(weights,rCtildeCalib,rBCalib);
            std::cout<<"RBI(0) loaded:"<<std::endl<<RBI<<std::endl;
            hasRBI=true;

            //also set initial reference position for CF
            rRefImu = internal_rI;
        }
    }

    double ttime=msg->tSolution.secondsOfWeek + msg->tSolution.fractionOfSecond + msg->tSolution.week * sec_in_week
          - msg->deltRSec;

    //if everything is working
    if(ttime>lastA2Dtime)  //Only use newest time. Ignore 0 messages.
    {
        hasAlreadyReceivedA2D=true;
        lastA2Dtime=ttime;
        if(msg->testStat > minTestStat)
        {
            validA2Dtest=true;
            //Store constrained baseline vector in I frame
            Eigen::Vector3d constrainedBaselineECEF(msg->rx,msg->ry,msg->rz);
            internal_rC = R_G2wrw*constrainedBaselineECEF;
            //internal_rC = -1*internal_rC; //match 1->2 sign convention if pprx reports 1->2
            internal_rC=unit3(internal_rC);
        }else
        {
            validA2Dtest=false;
        }

        //If the message is new, both messages have been received, and all teststats are good, then publish.
        if(abs(lastRTKtime-lastA2Dtime)<.001 && validA2Dtest && validRTKtest
                && hasAlreadyReceivedRTK && hasAlreadyReceivedA2D)  //only resend pose if new
        {
            internalSeq++;

            //Run CF
            Fimu = getFmatrixCF(ttime-tLastProcessed,imuAccelMeas,imuAttRateMeas,RBI) * Fimu;
            runCF(ttime-tLastProcessed);

            //Publish messages
            publishOdomAndMocap();

            //Clean up
            Fimu = Eigen::MatrixXd::Identity(15,15);
            tLastProcessed = ttime;
            RBI = updateRBIfromGamma(RBI,xState.middleRows(6,3));
            xState.middleRows(6,3)=Eigen::Vector3d::Zero();
            rRefImu=rRefImu+xState.topRows(3);
            xState.topRows(3)=Eigen::Vector3d::Zero();

            //Reset to avoid publishing twice
            hasAlreadyReceivedRTK=false; hasAlreadyReceivedA2D=false;
        }else //If not, tag this as correctly timed and pass off to RTK
        {
            lastA2Dtime=msg->tSolution.secondsOfWeek+msg->tSolution.fractionOfSecond-msg->deltRSec +
                        msg->tSolution.week * sec_in_week;
        }
    }
}


//Publish local_odom and mavros mocap
void gpsImuNode::publishOdomAndMocap()
{
    //not fully populated
    //std::cout << "publisher function called" << std::endl;
    nav_msgs::Odometry localOdom_msg;
    localOdom_msg.pose.pose.position.x = xState(0);
    localOdom_msg.pose.pose.position.y = xState(1);
    localOdom_msg.pose.pose.position.z = xState(2);
    localOdom_msg.twist.twist.linear.x = xState(3);
    localOdom_msg.twist.twist.linear.y = xState(4);
    localOdom_msg.twist.twist.linear.z = xState(5);
    Eigen::Vector3d gamma0(xState(6),xState(7),xState(8));
    Eigen::Quaterniond q0 = rotmat2quat(updateRBIfromGamma(RBI,gamma0));
    localOdom_msg.pose.pose.orientation.x=q0.x();
    localOdom_msg.pose.pose.orientation.y=q0.y();
    localOdom_msg.pose.pose.orientation.z=q0.z();
    localOdom_msg.pose.pose.orientation.w=q0.w();
    //Remove biases
    localOdom_msg.twist.twist.angular.x=imuAttRateMeas(0)-xState(12);
    localOdom_msg.twist.twist.angular.y=imuAttRateMeas(1)-xState(13);
    localOdom_msg.twist.twist.angular.z=imuAttRateMeas(2)-xState(14);
    //Publish local odometry message
    localOdom_pub_.publish(localOdom_msg);

    //px4 mocap topic to align frames
    geometry_msgs::PoseStamped mocap_msg;
    mocap_msg.pose.position = localOdom_msg.pose.pose.position;
    mocap_msg.pose.orientation = localOdom_msg.pose.pose.orientation;
//    mocap_msg.header = msg->header;
    mocap_msg.header.frame_id = "fcu";
    mocap_pub_.publish(mocap_msg);
}


//returns RRT. Must be converted to ORT (via tOffset) and then to gpsTime (via dtrx/clight)
void gpsImuNode::updateIMUtimeRRT(const uint64_t tIndex0, int &gpsWeek, int &gpsSec, float &gpsFracSec)
{

  static const int SF_TL = 24;
  static const int32_t SF_T = 0x1 << SF_TL;
  static const int32_t SF_T_MASK = SF_T - 1;
  static const int SEC_PER_WEEK(604800);
  static const double EPSILON(1e-15);
  
  //setWithTruncatedSampleTime()
  //float tFracIndex = tIndexKconfig;
  float tFracIndex = 0;
  //Modified slightly from gss->basetime.cpp
  //NOTE: truncL <=> msg->tIndexTrunc;
  //constexpr uint64_t one = 1ul; //defined at class level
  const uint64_t trunc = one << 39-14;  //=39-SF_SL
  const uint64_t truncHalf = trunc>>1;
  const uint64_t truncM = trunc-1;
  //trefWeek/FracSecs/SecOfWeek taken from NavigationSolution message
  //trefWeek=0; trefFracSecs=0; trefSecOfWeek=0;
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
  uint64_t fulltIndex = (reftIndex & ~truncM) | tIndex0;
  uint64_t truncReftIndex = reftIndex & truncM;
  if(truncReftIndex >= truncHalf){
    if(tIndex0 < truncReftIndex-truncHalf)
      fulltIndex += trunc;
  }
  else{
    if(tIndex0 > truncReftIndex+truncHalf)
      fulltIndex -= trunc;
  }

  //setWithSampleTime()
  uint64_t tIndex = fulltIndex;
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
  }

  /*std::cout << "Week: " << week_ << std::endl;
  std::cout << "Sec : " << secondsOfWeek_ << std::endl;
  std::cout << "Fsec: " << fractionOfSecond_ << std::endl;*/

  //return via non-const vars
  gpsWeek = week_;
  gpsSec = secondsOfWeek_;
  gpsFracSec = fractionOfSecond_;
  return;
}

} //end namespace
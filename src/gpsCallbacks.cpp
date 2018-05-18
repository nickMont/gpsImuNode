#include <Eigen/Geometry>
#include "estimationNode.hpp"
#include <string>
#include <iostream>
//Contains the two "wordy" A2D and SBRTK callbacks.
//Also contains some rotation matrix construction

namespace gpsimu_odom
{
void estimationNode::singleBaselineRTKCallback(const gbx_ros_bridge_msgs::SingleBaselineRTK::ConstPtr &msg)
{
  double ttime=msg->tSolution.secondsOfWeek + msg->tSolution.fractionOfSecond
      + msg->tSolution.week * sec_in_week - msg->deltRSec;
  if(ttime>lastRTKtime)  //only use newest time
  {
    hasAlreadyReceivedRTK=true;
    lastRTKtime=ttime;
    //If the message is accepted
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
      }else
      {
        validRTKtest=false;
      }

      //If the time is new, both messages for this time have been received, and teststats are good
    if(abs(lastRTKtime-lastA2Dtime)<.001 && validA2Dtest && validRTKtest
         && hasAlreadyReceivedRTK && hasAlreadyReceivedA2D)  //only resend pose if new
    {
      internalSeq++;
      double dtLastProc = ttime-tLastProcessed;
      //if(dtLastProc<0)
      //  { std::cout << "Negative times!" << "gps: " << ttime << "  imu: "<< tLastProcessed <<std::endl;}
      if(isCalibrated && dtLastProc>0)
      {
        runUKF(dtLastProc);

        //Publish messages
        P_report = Pimu;      
        updateType = "gps";
        publishOdomAndMocap();
  
        //Clean up
        Fimu = Eigen::MatrixXd::Identity(15,15);
        tLastProcessed = ttime;
      }else if(ttime-tProcPrev > 0 && isCalibrated && tryToUseBuffer) //If gps message is received out of order
      {
        runUKF_fromBuffer(ttime-tProcPrev);
        spkfPropagate15(xState,Pimu,Qk12dividedByDt*(tLastProcessed-ttime),tLastProcessed-ttime,
            imuAccelPrev,imuAttRatePrev,RBI,l_imu, Pimu,xState);
        //Do NOT publish
      }

      //Reset to avoid publishing twice
      hasAlreadyReceivedRTK=false; hasAlreadyReceivedA2D=false; 
    }
  }
}


void estimationNode::attitude2DCallback(const gbx_ros_bridge_msgs::Attitude2D::ConstPtr &msg)
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
        }
      return;
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
            const Eigen::Vector3d constrainedBaselineECEF(msg->rx,msg->ry,msg->rz);
            internal_rC = R_G2wrw*constrainedBaselineECEF;
            internal_rC = unit3(internal_rC);
        }else
        {
            validA2Dtest=false;
        }

        //If the message is new, both messages have been received, and all teststats are good, then publish.
        if(abs(lastRTKtime-lastA2Dtime)<.001 && validA2Dtest && validRTKtest
                && hasAlreadyReceivedRTK && hasAlreadyReceivedA2D)  //only resend pose if new
        {
          internalSeq++;
          double dtLastProc = ttime-tLastProcessed;
          //if(dtLastProc<0)
          //  { std::cout << "Negative times!" << "gps: " << ttime << "  imu: "<< tLastProcessed <<std::endl;}
          if(isCalibrated && dtLastProc>0)
          {
            runUKF(dtLastProc);

            //Publish messages
            P_report = Pimu;      
            updateType = "gps";
            publishOdomAndMocap();

            //Clean up
            Fimu = Eigen::MatrixXd::Identity(15,15);
            tLastProcessed = ttime;
          }else if(ttime-tProcPrev > 0 && isCalibrated && tryToUseBuffer) //If gps message is received out of order
          {
            runUKF_fromBuffer(ttime-tProcPrev);
            spkfPropagate15(xState,Pimu,Qk12dividedByDt*(tLastProcessed-ttime),tLastProcessed-ttime,
                imuAccelMeas,imuAttRateMeas,RBI,l_imu, Pimu,xState);
            //Do NOT publish
          }

          //Reset to avoid publishing twice
          hasAlreadyReceivedRTK=false; hasAlreadyReceivedA2D=false;
        }
    }
}


void estimationNode::navsolCallback(const gbx_ros_bridge_msgs::NavigationSolution::ConstPtr &msg)
{
  dtRX_meters = msg->deltatRxMeters;
}


//Get reference RRT time and measurement offset time from Observables message
void estimationNode::tOffsetCallback(const gbx_ros_bridge_msgs::ObservablesMeasurementTime::ConstPtr &msg)
{
  if(msg->tOffset.week<1e-9)
    {return;}
  toffsetWeek = msg->tOffset.week;
  toffsetSecOfWeek = msg->tOffset.secondsOfWeek;
  toffsetFracSecs = msg->tOffset.fractionOfSecond;
}


//Get upper 32 bits of tIndex counter
void estimationNode::imuConfigCallback(const gbx_ros_bridge_msgs::ImuConfig::ConstPtr &msg)
{
  ROS_INFO("Config message received.");
  imuConfigAccel = msg->lsbToMetersPerSecSq; //scaling to m/s2 from "non-engineering units"
  imuConfigAttRate = msg->lsbToRadPerSec; //scaling to rad/s from "non-engineering units"
  //imuSampleFreq = msg->sampleFreqNumerator/msg->sampleFreqDenominator/36/3600;  //samples per second
  
  sampleFreqNum = msg->sampleFreqNumerator;
  sampleFreqDen = msg->sampleFreqDenominator;
  tIndexConfig = msg->tIndexk;
}


//Publish local_odom and mavros mocap
void estimationNode::publishOdomAndMocap()
{

    //Update rotation matrix and force orthonormality  
    RBI=updateRBIfromGamma(RBI, xState.middleRows(6,3));
    xState.middleRows(6,3)=Eigen::Vector3d::Zero();
    //RBI = orthonormalize(RBI);

    /*if(!updateType.compare("gps"))
      {std::cout << "UPDATE TYPE: "<< updateType<< " " << updateType<< " " << updateType<< " " << updateType<< " " << updateType <<std::endl;}
    else{
    std::cout << "UPDATE TYPE: " << updateType <<std::endl;}*/
    //std::cout << "Accelbias:" <<std::endl;
    //std::cout << xState(9) << " " << xState(10) << " " << xState(11) << std::endl;
    //std::cout << "P_eigs: " <<std::endl << (Pimu.eigenvalues()).transpose() << std::endl;

    //std::cout << "Gyrobias x100:" <<std::endl;
    //std::cout << xState(12)*100.0 << " " << xState(13)*100.0 << " " << xState(14)*100.0 << std::endl;

    //Output counter
    static int counter(0);
    counter++;
    if(counter%10==0)
    {
      RBI = orthonormalize(RBI);
    }
    //{std::cout<<"RBI:"<<std::endl<<RBI<<std::endl;}
    
    //not fully populated
    //std::cout << "publisher function called" << std::endl;
    nav_msgs::Odometry localOdom_msg;
    //std::cout << "xk:" << std::endl << xState.topRows(3) << std::endl;

    //Eigen::Matrix<double,6,6> covXG;
    //covXG.topLeftCorner(3,3) = P_report.topLeftCorner(3,3);
    //covXG.bottomRightCorner(3,3) = P_report.middleRows(6,6,3,3);
    //covXG.topRightCorner(3,3) = P_report.block(0,6,3,3);
    //covXG.bottomLeftCorner(3,3) = P_report.block(6,0,3,3);

    //Eigen::Matrix3d H;
    //H<<1,?,?

    //Generate message
    localOdom_msg.header.stamp = ros::Time::now();
    localOdom_msg.header.frame_id = updateType;
    localOdom_msg.child_frame_id = "world";
    localOdom_msg.pose.pose.position.x = xState(0);
    localOdom_msg.pose.pose.position.y = xState(1);
    localOdom_msg.pose.pose.position.z = xState(2);
    localOdom_msg.twist.twist.linear.x = xState(3);
    localOdom_msg.twist.twist.linear.y = xState(4);
    localOdom_msg.twist.twist.linear.z = xState(5);
    Eigen::Quaterniond q0 = rotmat2quat(RBI);
    localOdom_msg.pose.pose.orientation.x=q0.x();
    localOdom_msg.pose.pose.orientation.y=q0.y();
    localOdom_msg.pose.pose.orientation.z=q0.z();
    localOdom_msg.pose.pose.orientation.w=q0.w();
    //Remove biases
    localOdom_msg.twist.twist.angular.x=imuAttRateMeas(0)-xState(12);
    localOdom_msg.twist.twist.angular.y=imuAttRateMeas(1)-xState(13);
    localOdom_msg.twist.twist.angular.z=imuAttRateMeas(2)-xState(14);
    /*for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
      {
        localOdom_msg.pose.covariance[6*i + j] = P_report(i, j);
        localOdom_msg.twist.covariance[6*i + j] = P_report(3+i, 3+j);
        //Covariance of attitude rate is gaussian(IMU error) + gaussian(bias estimate)
        localOdom_msg.twist.covariance[6*i + j + 21] = Qimu(3+i, 3+j) + P_report(12+i, 12+j);
      }
    }*/

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


} //end namespace
#include <Eigen/Geometry>
#include "estimationNode.hpp"
#include <string>
#include <iostream>
//Contains the two "wordy" A2D and SBRTK callbacks.
//Also contains some rotation matrix construction

namespace gpsimu_odom
{
    
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
    const long long int tIndexFull = (tIndexConfig_ & ~mask) | (tIndex & mask); //Bit-mask, don't bit-shift
    const double sampleFreq = sampleFreqNum_/sampleFreqDen_;
    const double tRRT = tIndexFull/sampleFreq; //tIndexFull is in samples, divide by samples/s
    const double tOffset = toffsetFracSecs_ + toffsetWeek_*SEC_PER_WEEK + toffsetSecOfWeek_;
    const double tORT = tRRT + tOffset; //tOffset comes from ObservablesMeasurementTime
    const double tGPS = tORT - dtRXinMeters_/cLight; //dtrx comes from NavigationSolution
    const double thisTime = tGPS;
    std::cout.precision(17);
    dt = thisTime-tLastImu;
    if(dt<=1e-9) //Note: NOT abs(dt), as this checks whether or not messages are received out of order as well
    {
        std::cout << "Error: 1ns between IMU measurements" << std::endl;
        //std::cout << "Time: " << thisTime-double(toffsetWeek_*SEC_PER_WEEK) <<std::endl;
        return;
    }
    //Only update last time used IF this time is accepted
    tLastImu=thisTime;

    //Use scale parameter
    Eigen::Vector3d imuAccelMeas, imuAttRateMeas;
    imuAccelMeas(0) = msg->acceleration[0] * imuConfigAccel_;
    imuAccelMeas(1) = msg->acceleration[1] * imuConfigAccel_;
    imuAccelMeas(2) = msg->acceleration[2] * imuConfigAccel_;
    imuAttRateMeas(0) = msg->angularRate[0] * imuConfigAttRate_;
    imuAttRateMeas(1) = msg->angularRate[1] * imuConfigAttRate_;
    imuAttRateMeas(2) = msg->angularRate[2] * imuConfigAttRate_;

    //Rotate gyro/accel measurements to body frame
    //NOTE1: this does not need Rpqr convention
    //NOTE2: this is hardcoded for the lynx as mounted on phoenix and company
    Eigen::Matrix3d Raccel, Rgyro;
    Raccel << 0,-1,0, -1,0,0, 0,0,-1;
    Rgyro  << 0,-1,0, -1,0,0, 0,0,-1;
    imuAccelMeas = Raccel*imuAccelMeas;
    imuAttRateMeas = Rgyro*imuAttRateMeas;

    //Run CF if calibrated
    if(isCalibrated_)
    {
        imuSeq++;
        double dtLastProc = thisTime - tLastProcessed_;
        if(dtLastProc>0) //Just in case a gps message is received late
        { 
            //Construct measurements
            const imuMeas thisImuMeas(thisTime,imuAccelMeas,imuAttRateMeas);
            lastImuMeas_ = thisImuMeas;

            imuFilter_.runUKFpropagateOnly(tLastProcessed_,thisImuMeas);

            //Publish
            updateType = "imu"; publishOdomAndMocap();
            //RBI_ is updated when publisher is called.

            //Cleanup
            tLastProcessed_ = thisTime;

            //Warn if signal is lost
            if( (thisTime-lastRTKtime_ > 0.5) || (thisTime-lastA2Dtime_ > 0.5) )
            {
                if(warnCounter%40==0)
                {ROS_INFO("GPS outage warning!");}
                warnCounter++;
            }else
            {
                if(warnCounter!=0)
                {ROS_INFO("GPS restored.");}
                warnCounter=0;
            }
        }
    }else if(rbiIsInitialized_)  //if RBI has been calculated but the biases have not been calculated
    { 
        ba0=ba0+0.01*(imuAccelMeas - RBI_.transpose()*Eigen::Vector3d(0,0,9.8)); //inefficient
        bg0=bg0+0.01*imuAttRateMeas;
        Eigen::Vector3d rI0;
        rI0 = rPrimaryMeas_ - RBI_.transpose()*Lcg2p_;
        Eigen::Matrix<double,15,1> xState;
        xState<<rI0(0),rI0(1),rI0(2), 0,0,0, 0,0,0, ba0(0),ba0(1),ba0(2), bg0(0),bg0(1),bg0(2);
        counter++;
        // Try ground calibration step for simplicity
        if(counter>=100)
        {
            isCalibrated_ = true;
            tLastProcessed_ = thisTime;
            imuFilter_.setState(xState,RBI_);
            //std::cout << "state0" <<std::endl<<xState<<std::endl;
        }
    }

}

void estimationNode::singleBaselineRTKCallback(const gbx_ros_bridge_msgs::SingleBaselineRTK::ConstPtr &msg)
{
    double ttime=msg->tSolution.secondsOfWeek + msg->tSolution.fractionOfSecond
            + msg->tSolution.week * sec_in_week_ - msg->deltRSec;
    if(ttime>lastRTKtime_)  //only use newest time
    {
        hasAlreadyReceivedRTK_=true;
        lastRTKtime_=ttime;
        //If the message is accepted
        if(msg->testStat > minTestStat_)
            {
                validRTKtest_=true;
                Eigen::Vector3d tmpvec;
                //Rotate rECEF to rI and store in rPrimaryMeas_
                tmpvec(0) = msg->rxRov - zeroInECEF_(0); //error vector from ECEF at init time
                tmpvec(1) = msg->ryRov - zeroInECEF_(1);
                tmpvec(2) = msg->rzRov - zeroInECEF_(2);
                rPrimaryMeas_ = Rwrw_*(Recef2enu_*tmpvec - zeroInWRW_);
                //ROS_INFO("%f %f %f %f",msg->rx, msg->rxRov, tmpvec(0), internalPose(0)); //debugging
            }else
            {
                validRTKtest_=false;
            }

            //If the time is new, both messages for this time have been received, and teststats are good
        if(abs(lastRTKtime_-lastA2Dtime_)<.001 && validA2Dtest_ && validRTKtest_
                 && hasAlreadyReceivedRTK_ && hasAlreadyReceivedA2D_)  //only resend pose if new
        {
            internalSeq++;
            double dtLastProc = ttime-tLastProcessed_;

            if(isCalibrated_ && dtLastProc>0)
            {
                const gpsMeas thisGpsMeas(ttime,rPrimaryMeas_,rS2PMeas_);
                imuFilter_.runUKF(lastImuMeas_,thisGpsMeas);

                //Publish messages
                updateType = "gps";
                publishOdomAndMocap();
                //RBI_ is updated when publisher is called.

                tLastProcessed_ = ttime;
            }

            //Reset to avoid publishing twice
            hasAlreadyReceivedRTK_=false; hasAlreadyReceivedA2D_=false; 
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

    if(!rbiIsInitialized_ && msg->testStat>=100)
    {
        const Eigen::Vector3d constrainedBaselineECEF(msg->rx,msg->ry,msg->rz);
        const Eigen::Vector3d constrainedBaselineI(unit3(Rwrw_*Recef2enu_*constrainedBaselineECEF));
        rCtildeCalib_(rCCalibCounter%calibSamples,0)=constrainedBaselineI(0);
        rCtildeCalib_(rCCalibCounter%calibSamples,1)=constrainedBaselineI(1);
        rCtildeCalib_(rCCalibCounter%calibSamples,2)=constrainedBaselineI(2);
        rBCalib_(rCCalibCounter%calibSamples,0)=1;
        rBCalib_(rCCalibCounter%calibSamples,1)=0;
        rBCalib_(rCCalibCounter%calibSamples,2)=0;
        //rCB=???
        rCCalibCounter++;

        //estimate initial RBI from wahba
        if(rCCalibCounter>=calibSamples)
        {
            rCtildeCalib_(calibSamples,0)=0; rBCalib_(calibSamples,0)=0;
            rCtildeCalib_(calibSamples,1)=0; rBCalib_(calibSamples,1)=0;
            rCtildeCalib_(calibSamples,2)=1; rBCalib_(calibSamples,2)=1;
            Eigen::MatrixXd weights;
            weights.resize(calibSamples+1,1);
            weights.topRows(calibSamples)=0.5*1/calibSamples*Eigen::MatrixXd::Ones(calibSamples,1);
            weights(calibSamples)=0.5;
            RBI_=rotMatFromWahba(weights,rCtildeCalib_,rBCalib_);
            rbiIsInitialized_=true;
        }
        return;
    }


    double ttime=msg->tSolution.secondsOfWeek + msg->tSolution.fractionOfSecond + msg->tSolution.week * sec_in_week_
        - msg->deltRSec;

    //if everything is working
    if(ttime>lastA2Dtime_)  //Only use newest time. Ignore 0 messages.
    {
        hasAlreadyReceivedA2D_=true;
        lastA2Dtime_=ttime;
        if(msg->testStat > minTestStat_)
        {
            validA2Dtest_=true;
            //Store constrained baseline vector in I frame
            const Eigen::Vector3d constrainedBaselineECEF(msg->rx,msg->ry,msg->rz);
            rS2PMeas_ = Recef2wrw_*constrainedBaselineECEF;
            rS2PMeas_ = unit3(rS2PMeas_);
        }else
        {
            validA2Dtest_=false;
    }

    //If the message is new, both messages have been received, and all teststats are good, then publish.
    if(abs(lastRTKtime_-lastA2Dtime_)<.001 && validA2Dtest_ && validRTKtest_
        && hasAlreadyReceivedRTK_ && hasAlreadyReceivedA2D_)  //only resend pose if new
    {
        internalSeq++;
        double dtLastProc = ttime-tLastProcessed_;
        //if(dtLastProc<0)
        //  { std::cout << "Negative times!" << "gps: " << ttime << "  imu: "<< tLastProcessed_ <<std::endl;}
        if(isCalibrated_ && dtLastProc>0)
        {
            const gpsMeas thisGpsMeas(ttime,rPrimaryMeas_,rS2PMeas_);
            imuFilter_.runUKF(lastImuMeas_,thisGpsMeas);

            //Publish messages
            updateType = "gps";
            publishOdomAndMocap();
            //RBI_ is updated when publisher is called.

            tLastProcessed_ = ttime;
        }

        //Reset to avoid publishing twice
        hasAlreadyReceivedRTK_=false; hasAlreadyReceivedA2D_=false;
        }
    }
}


void estimationNode::navsolCallback(const gbx_ros_bridge_msgs::NavigationSolution::ConstPtr &msg)
{
    dtRXinMeters_ = msg->deltatRxMeters;
}


//Get reference RRT time and measurement offset time from Observables message
void estimationNode::tOffsetCallback(const gbx_ros_bridge_msgs::ObservablesMeasurementTime::ConstPtr &msg)
{
    if(msg->tOffset.week<1e-9)
        {return;}
    toffsetWeek_ = msg->tOffset.week;
    toffsetSecOfWeek_ = msg->tOffset.secondsOfWeek;
    toffsetFracSecs_ = msg->tOffset.fractionOfSecond;
}


//Get upper 32 bits of tIndex counter
void estimationNode::imuConfigCallback(const gbx_ros_bridge_msgs::ImuConfig::ConstPtr &msg)
{
    ROS_INFO("Config message received.");
    imuConfigAccel_ = msg->lsbToMetersPerSecSq; //scaling to m/s2 from "non-engineering units"
    imuConfigAttRate_ = msg->lsbToRadPerSec; //scaling to rad/s from "non-engineering units"
    //imuSampleFreq = msg->sampleFreqNumerator/msg->sampleFreqDenominator/36/3600;  //samples per second
    
    sampleFreqNum_ = msg->sampleFreqNumerator;
    sampleFreqDen_ = msg->sampleFreqDenominator;
    tIndexConfig_ = msg->tIndexk;
}


//Publish local_odom and mavros mocap
void estimationNode::publishOdomAndMocap()
{
    Eigen::Matrix<double,15,15> Preport;
    imuFilter_.getCovariance(Preport);

    //Update rotation matrix and force orthonormality  
    Eigen::Matrix<double,15,1> xState;
    Eigen::Vector3d imuAccel,imuAttRate;
    double tlastImu;
    imuFilter_.getState(xState, RBI_);

    nav_msgs::Odometry localOdom_msg;

    //Generate message
    localOdom_msg.header.stamp = ros::Time::now();
    localOdom_msg.header.frame_id = "world";
    localOdom_msg.child_frame_id = "world";
    localOdom_msg.pose.pose.position.x = xState(0);
    localOdom_msg.pose.pose.position.y = xState(1);
    localOdom_msg.pose.pose.position.z = xState(2);
    localOdom_msg.twist.twist.linear.x = xState(3);
    localOdom_msg.twist.twist.linear.y = xState(4);
    localOdom_msg.twist.twist.linear.z = xState(5);
    //Eigen::Quaterniond q0 = rotmat2quat(RBI_);
    Eigen::Quaterniond q0(RBI_);
    localOdom_msg.pose.pose.orientation.x=q0.x();
    localOdom_msg.pose.pose.orientation.y=q0.y();
    localOdom_msg.pose.pose.orientation.z=q0.z();
    localOdom_msg.pose.pose.orientation.w=q0.w();

    //Remove biases
    lastImuMeas_.getMeas(tlastImu,imuAccel,imuAttRate);
    localOdom_msg.twist.twist.angular.x=imuAttRate(0)-xState(12);
    localOdom_msg.twist.twist.angular.y=imuAttRate(1)-xState(13);
    localOdom_msg.twist.twist.angular.z=imuAttRate(2)-xState(14);
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            localOdom_msg.pose.covariance[6*i + j] = Preport(i, j);
            localOdom_msg.twist.covariance[6*i + j] = Preport(3+i, 3+j);
            //Covariance of attitude rate is gaussian(IMU error) + gaussian(bias estimate)
            localOdom_msg.twist.covariance[6*i + j + 21] = QgyroOutput_(i, j) + Preport(12+i, 12+j);
        }
    }

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
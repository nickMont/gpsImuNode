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
#include <stdint.h>
#include <cmath>

//Contains data storage classes

class gpsTime
{
  	public:
    	gpsTime(){week_=0.0;sec_=0.0;fracSec_=0.0;}
    	gpsTime(const double wk, const double se, const double fSec){week_=wk;sec_=se;fracSec_=fSec;}
    	void getTime(double &wk, double &se, double &fSec){wk=week_;se=sec_;fSec=fracSec_;}
    	void setTime(const double wk, const double se, const double fSec){week_=wk;sec_=se;fracSec_=fSec;}
    	double toSec(){return (week_*604800.0+sec_+fracSec_);}
  	private:
    	double week_, sec_, fracSec_;
};


//All measurements are expressed in ECEF
class gpsMeas
{
	public:
		gpsMeas() {tSec_=0.0; rPrimary_=Eigen::Vector3d::Zero(); rS2P_=Eigen::Vector3d::Zero();
			Recef2wrw_=Eigen::Matrix3d::Identity();}
		gpsMeas(const double t, const Eigen::Vector3d &rP, const Eigen::Vector3d &rC, const Eigen::Matrix3d &R)
			{tSec_=t; rPrimary_=rP; rS2P_=rC; Recef2wrw_=R;};
		void setRotMats(const Eigen::Matrix3d &Rwrw, const Eigen::Matrix3d &Recef2enu)
			{Recef2wrw_ = Rwrw * Recef2enu; }
		void setMeas(const double t, const Eigen::Vector3d &rP, const Eigen::Vector3d &rC)
			{tSec_=t; rPrimary_=rP; rS2P_=rC;}
		void getMeasEcef(double &t, Eigen::Vector3d &rP, Eigen::Vector3d &rC)
			const {t=tSec_; rP=rPrimary_; rC=rS2P_;}
		void getMeasEnu(double &t, Eigen::Vector3d &rI, Eigen::Vector3d &rC)
			const {t=tSec_; rI=Recef2wrw_*rPrimary_; rC=Recef2wrw_*rS2P_;}
		void getTime(double &t) const {t=tSec_;}
	private:
		double tSec_;
		Eigen::Vector3d rPrimary_, rS2P_;
		Eigen::Matrix3d Recef2wrw_;
		// vWRW = Recef2wrw_ * vECEF
};

class imuMeas
{
	public:
		imuMeas(){tSec_=0.0; a_=Eigen::Vector3d::Zero(); wB_=Eigen::Vector3d::Zero();}
		imuMeas(const double t, const Eigen::Vector3d &aa, const Eigen::Vector3d &ww)
			{tSec_=t; a_=aa; wB_=ww;};
		void setMeas(const double t, const Eigen::Vector3d &aIn, const Eigen::Vector3d &wIn)
			{tSec_=t; a_=aIn; wB_=wIn;};
		void getMeas(double &t, Eigen::Vector3d &aOut, Eigen::Vector3d &wOut)
			const {t=tSec_; aOut=a_; wOut=wB_;}
		void getTime(double &t) const {t=tSec_;}	
	private:
		double tSec_;
		Eigen::Vector3d a_, wB_;
};


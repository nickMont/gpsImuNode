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

#include "filter.h"
#include "filterTW.h"

namespace gpsimu_odom
{
class gpsImu
{
	public:
      	Eigen::Matrix3d updateRBIfromGamma(const Eigen::Matrix3d R0, const Eigen::Vector3d gamma);
      	Eigen::Matrix3d hatmat(const Eigen::Vector3d v1);
        Eigen::Matrix3d rotMatFromEuler(Eigen::Vector3d ee);
        Eigen::Matrix3d rotMatFromWahba(const Eigen::VectorXd weights, const::Eigen::MatrixXd vI,
    	   const::Eigen::MatrixXd vB);
        Eigen::Quaterniond rotmat2quat(const Eigen::Matrix3d RR);
        Eigen::Vector3d unit3(const Eigen::Vector3d v1);
        Eigen::Matrix3d orthonormalize(const Eigen::Matrix3d inmat);
        Eigen::Matrix3d rotMatFromQuat(const Eigen::Quaterniond qq);
        void saturateBiases(const double baMax, const double bgMax);
        double symmetricSaturationDouble(const double inval, const double maxval);
        Eigen::Matrix3d euler2dcm312(const Eigen::Vector3d ee);
        void spkfPropagate15(const Eigen::Matrix<double,15,1> *x0, const Eigen::Matrix<double,15,15> P0,
            const Eigen::Matrix<double,12,12> Q, const double dt, const Eigen::Vector3d fB0, const Eigen::Vector3d wB0,
            const Eigen::Matrix3d RR, const Eigen::Vector3d lAB, Eigen::Matrix<double,15,15> &Pkp1, Eigen::Matrix<double,15,1> &xkp1);
        void spkfMeasure6(const Eigen::Matrix<double,15,1> *x0, const Eigen::Matrix<double,15,15> P0,
            const Eigen::Matrix<double,6,6> R, const Eigen::Vector3d rI_measurement, const Eigen::Vector3d rCu_measurement,
    	    const Eigen::Matrix3d RR, const Eigen::Vector3d lcg2p, const Eigen::Vector3d ls2p,
    	    Eigen::Matrix<double,15,15> &Pkp1, Eigen::Matrix<double,15,1> &xkp1);
  	    Eigen::Matrix<double,15,15> getNumderivF(const double dv, const double dt,
    	    const Eigen::Matrix<double,15,1> x0,const Eigen::Vector3d fB0, const Eigen::Vector3d wB0,
    	    const Eigen::Matrix3d RR, const Eigen::Vector3d lAB);
  	    Eigen::Matrix<double,6,1> hnonlinSPKF(const Eigen::Matrix<double,15,1> *x0,
     	    const Eigen::Matrix3d RR, const Eigen::Vector3d ls2p, const Eigen::Vector3d lcg2p,
   	 	    const Eigen::Matrix<double,6,1> vk);
  	    Eigen::Matrix<double,15,1> fdynSPKF(const Eigen::Matrix<double,15,1> *x0, const double dt,
    	   const Eigen::Vector3d fB0, const Eigen::Matrix<double,12,1> vk, const Eigen::Vector3d wB0,
    	   const Eigen::Matrix3d RR, const Eigen::Vector3d lAB);
  	    void runUKF(double dt0); 
	    void runUKF_fromBuffer(double dt0);


	private:

	 Eigen::Vector3d baseECEF_vector, baseENU_vector, WRW0_ecef, arenaRefCenter,
    	internal_rI, internal_rC, internal_rImu, rRefImu, n_err, imuAccelMeas,
    	imuAttRateMeas, l_imu, l_s2p, l_cg2p;

  	Eigen::Matrix3d Recef2enu, Rwrw, R_G2wrw, RBI;
  	Eigen::Matrix<double,21,3> rCtildeCalib, rBCalib;
  	Eigen::Matrix<double,15,1> xState, xStatePrev;
  	//Pimu, Fimu are P, F matrices of EKF.  P_report is special IMU component used in mocap for publisher.
  	//Pimu sets P_report; P_report is used exclusively for reporting P.
  	Eigen::Matrix<double,15,15> Fimu, Pimu, P_report, PimuPrev;
  	Eigen::Matrix<double,6,6> Qimu, Rk;

  	Eigen::Quaterniond internalQuat, quaternionSetpoint;
  	int centerFlag, internalSeq;
  	double lastRTKtime, lastA2Dtime, minTestStat, max_accel, throttleSetpoint, throttleMax,
    	imuConfigAccel, imuConfigAttRate, tMeasOffset, pi, tLastProcessed;

    long long int tIndexConfig;
  	Eigen::Matrix<double,12,12> Qk12, Qk12dividedByDt;
  	uint32_t tmeasWeek, tmeasSecOfWeek, toffsetWeek, toffsetSecOfWeek;
  	double tmeasFracSecs, toffsetFracSecs, tauG, tauA;
  	uint64_t sampleFreqNum, sampleFreqDen;
  	uint64_t one, imuTimeTrunc;
  	double dtRX_meters;
  	Eigen::Vector3d attRateMeasOrig, accelMeasOrig, initBA, initBG, imuAccelPrev, imuAttRatePrev;
  	int32_t tnavsolWeek, tnavsolSecOfWeek;
	double dtGPS, tnavsolFracSecs, maxBa, maxBg, sec_in_week, tProcPrev;
};

}
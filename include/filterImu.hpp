#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>
#include "classes.h"
#include <ros/ros.h>

namespace gpsimu_odom
{
class gpsImu
{
	public:
        gpsImu() {ROS_INFO("Remember to intialize IMU components individually after initial calibration");
            xState_=Eigen::Matrix<double,15,1>::Zero(); Pimu_=Eigen::Matrix<double,15,15>::Zero();}

        //Helper functions
      	Eigen::Matrix3d updateRBIfromGamma(const Eigen::Matrix3d &R0, const Eigen::Vector3d &gamma);
      	Eigen::Matrix3d hatmat(const Eigen::Vector3d &v1);
        Eigen::Matrix3d rotMatFromEuler(const Eigen::Vector3d &ee);
        Eigen::Matrix3d rotMatFromWahba(const Eigen::VectorXd &weights, const::Eigen::MatrixXd &vI,
    	   const::Eigen::MatrixXd &vB);
        Eigen::Quaterniond rotmat2quat(const Eigen::Matrix3d &RR);
        Eigen::Vector3d unit3(const Eigen::Vector3d &v1);
        Eigen::Matrix3d orthonormalize(const Eigen::Matrix3d &inmat);
        Eigen::Matrix3d rotMatFromQuat(const Eigen::Quaterniond &qq);
        void saturateBiases(const double baMax, const double bgMax);
        double symmetricSaturationDouble(const double inval, const double maxval);
        Eigen::Matrix3d euler2dcm312(const Eigen::Vector3d &ee);
        
        //UKF functions
        void spkfPropagate15(const Eigen::Matrix<double,15,1> &x0, const Eigen::Matrix<double,15,15> &P0,
            const Eigen::Matrix<double,12,12> &Q, const double dt, const Eigen::Vector3d &fB0, const Eigen::Vector3d &wB0,
            const Eigen::Matrix3d &RR, const Eigen::Vector3d &lAB, Eigen::Matrix<double,15,15> &Pkp1, Eigen::Matrix<double,15,1> &xkp1);
        void spkfMeasure6(const Eigen::Matrix<double,15,1> &x0, const Eigen::Matrix<double,15,15> &P0,
            const Eigen::Matrix<double,6,6> &R, const Eigen::Vector3d &rI_measurement, const Eigen::Vector3d &rCu_measurement,
    	    const Eigen::Matrix3d &RR, const Eigen::Vector3d &lcg2p, const Eigen::Vector3d &ls2p,
    	    Eigen::Matrix<double,15,15> &Pkp1, Eigen::Matrix<double,15,1> &xkp1);
  	    Eigen::Matrix<double,6,1> hnonlinSPKF(const Eigen::Matrix<double,15,1> &x0,
     	    const Eigen::Matrix3d &RR, const Eigen::Vector3d &ls2p, const Eigen::Vector3d &lcg2p,
   	 	    const Eigen::Matrix<double,6,1> &vk);
  	    Eigen::Matrix<double,15,1> fdynSPKF(const Eigen::Matrix<double,15,1> &x0, const double dt,
    	   const Eigen::Vector3d &fB0, const Eigen::Matrix<double,12,1> &vk, const Eigen::Vector3d &wB0,
    	   const Eigen::Matrix3d &RR, const Eigen::Vector3d &lAB);
  	    void runUKF(const imuMeas &imu, const gpsMeas &gps); 
        void runUKFpropagateOnly(const double lastTime, const imuMeas &imu);

        //Getters and setters. Everything must be initialized piecemeal. This is preferred as some things are not known
        //at initialization time.
        void setLevers(const Eigen::Vector3d &Lcg2p, const Eigen::Vector3d &Lcg2imu, const Eigen::Vector3d &Ls2p)
            {Lcg2p_=Lcg2p; Limu_=Lcg2imu; Ls2p_=Ls2p;}
        void setCovariances(const Eigen::Matrix<double,12,12> &Qin, const Eigen::Matrix<double,15,15> &Pin,
            const Eigen::Matrix<double,6,6> &Rin) {Qk12_=Qin; Pimu_=Pin; Rk_=Rin;}
        void setCovariances(const double dt, const Eigen::Matrix<double,12,12> &Qin, const Eigen::Matrix<double,15,15> &Pin,
            const Eigen::Matrix<double,6,6> &Rin) {Qk12_=Qin; Pimu_=Pin; Rk_=Rin; Qk12dividedByDt_=Qin/dt;}
        void setImuParams(const double tauA_in, const double tauG_in)
            {tauG_=tauG_in; tauA_=tauA_in;}
        void setState(const Eigen::Matrix<double,15,1> xIn, const Eigen::Matrix3d &RR)
            {xState_=xIn; RBI_=RR;}
        void setBiasSaturationLimits(const double baIn, const double bgIn) {maxBa_=baIn; maxBg_=bgIn;}
        void getCovariance(Eigen::Matrix<double,15,15> &Pout) {Pout=Pimu_;}
        void getState(Eigen::Vector3d &pos, Eigen::Vector3d &vel, Eigen::Matrix3d &RR)
            {pos=xState_.topRows(3); vel=xState_.middleRows(3,3); RR=RBI_;}
        void getState(Eigen::Matrix<double,15,1> &state, Eigen::Matrix3d &RR)
            {state=xState_; RR=RBI_;}

	private:
	    Eigen::Vector3d Limu_, Ls2p_, Lcg2p_;

      	Eigen::Matrix3d Recef2enu, Rwrw, R_G2wrw, RBI_;
        Eigen::Matrix<double,21,3> rCtildeCalib, rBCalib;
        Eigen::Matrix<double,15,1> xState_;
        //Pimu_, Fimu are P, F matrices of EKF.  P_report is special IMU component used in mocap for publisher.
        //Pimu_ sets P_report; P_report is used exclusively for reporting P. 
        Eigen::Matrix<double,15,15> Pimu_, P_report, Pimu_Prev;
        Eigen::Matrix<double,6,6> Qimu, Rk_;

     	Eigen::Matrix<double,12,12> Qk12_, Qk12dividedByDt_;
      	double tauG_, tauA_, maxBa_, maxBg_;
      	uint64_t one;
};

}
#include "gpsImuNode.hpp"
namespace gpsimu_odom
{

Eigen::Matrix<double,15,15> gpsImuNode::getFmatrixCF(const double dt, const Eigen::Vector3d fB,
	const Eigen::Vector3d omegaB, const Eigen::Matrix3d RBI)
{
	//A is continuous time, Fk is discrete time
	Eigen::Matrix<double,15,15> Ak = Eigen::MatrixXd::Zero(15,15);
	Ak.block(3,0,3,3)=Eigen::Matrix3d::Identity();
	Ak.block(3,6,3,3)=-RBI*hatmat(fB);
	Ak.block(3,12,3,3)=RBI;
	Ak.block(6,6,3,3)=hatmat(omegaB);
	return Eigen::Matrix<double,15,15>::Identity() + dt*Ak;
}

Eigen::Matrix<double,15,6> gpsImuNode::getGammakmatrixCF(const double dt, const Eigen::Matrix3d RBI)
{
	Eigen::Matrix<double,15,6> Gammak = Eigen::Matrix<double,15,6>::Zero();
	Eigen::Matrix3d eye3 = Eigen::Matrix3d::Identity();
	Gammak.block(0,3,3,3)=dt*eye3;
	Gammak.block(3,0,3,3)=dt*dt*0.5*RBI;
	Gammak.block(6,0,3,3)=dt*RBI;
	Gammak.block(9,3,3,3)=eye3;
	Gammak.block(12,0,3,3)=eye3;
	return Gammak;
}

Eigen::Matrix<double,3,15> gpsImuNode::getHkmatrixCF(const Eigen::Vector3d Lab, const Eigen::Matrix3d RBI)
{
	Eigen::Matrix<double,3,15> Hk = Eigen::Matrix<double,3,15>::Zero();
	Hk.block(0,0,3,3)=Eigen::Matrix3d::Identity();
	Hk.block(6,0,3,3)=-RBI*hatmat(Lab);
	return Hk;
}

//Cross product equivalent.  Named this way for consistency with nn_imu_dat
Eigen::Matrix3d gpsImuNode::hatmat(const Eigen::Vector3d v1)
{
	Eigen::Matrix3d f = Eigen::MatrixXd::Zero(3,3);
	f(0,1)=-v1(2); f(0,2)=v1(1);
	f(1,0)=v1(2); f(1,2)=-v1(0);
	f(2,0)=-v1(1); f(2,1)=v1(0);
	return f;
}






}

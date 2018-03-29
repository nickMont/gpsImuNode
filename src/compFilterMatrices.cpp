#include "gpsImuNode.hpp"
#include <Eigen/SVD>


namespace gpsimu_odom
{
//Runs the CF
void gpsImuNode::runCF(double dt0)
{
	std::cout<<"RUNNING CF, dt="<<dt0<<std::endl;
	Eigen::Matrix<double,15,15> pbar0;
	Eigen::Matrix<double,15,1>  xbar0;
	//construct measurement deltas
	//error state filter
	//internal_rImu = xState.topRows(3)+rRefImu;
	//Eigen::Vector3d drI = internal_rImu - internal_rI;
	//true state filter
	Eigen::Vector3d drI = internal_rI;
	//run EKF
	kfCFPropagate(dt0,Pimu,xState,Fimu,RBI,Qimu, pbar0,xbar0);
	kfCFMeasure2(xbar0, drI, internal_rC, RBI, pbar0, l_cg2p, l_s2p, Rk, Pimu, xState);
	//std::cout << "xstate after measurement" << std::endl << xState <<std::endl;
}


//CF propagation step
void gpsImuNode::kfCFPropagate(const double dt0, const Eigen::Matrix<double,15,15> P0,
    const Eigen::Matrix<double,15,1> x0, const Eigen::Matrix<double,15,15> F0,
    const Eigen::Matrix3d RR, const Eigen::Matrix<double,6,6> Qk,
    Eigen::Matrix<double,15,15> &Pbar, Eigen::Matrix<double,15,1> &xBar)
{
	Eigen::Matrix<double,15,6> gammak=getGammakmatrixCF(dtGPS,RR); //dt for gammak is from last gps to current gps
	xBar = fdyn(x0,dt0,imuAccelMeas,imuAttRateMeas,RR,l_imu);
	Pbar = F0*P0*F0.transpose()+gammak*Qk*gammak.transpose();
}


//CF measurement step for two antennas
void gpsImuNode::kfCFMeasure2(const Eigen::Matrix<double,15,1> xBar, const Eigen::Vector3d drI, const Eigen::Vector3d drS2P,
    const Eigen::Matrix3d RR, const Eigen::Matrix<double,15,15> Pbar,
    const Eigen::Vector3d Lcg2p, const Eigen::Vector3d Ls2p, const Eigen::Matrix<double,6,6> Rbar,
    Eigen::Matrix<double,15,15> &Pkp1, Eigen::Matrix<double,15,1> &xkp1)
{
	Eigen::Matrix<double,6,1> zk;
	zk.topRows(3)=drI;
	zk.bottomRows(3)=drS2P;
	//error state
	//Eigen::Matrix<double,6,15> Hk = getHkmatrixTwoAntennaCF(Limu,Ls2p,RR);
	//true state
	//std::cout << "xbar0: " << std::endl<<xBar<<std::endl;
	//std::cout << "zk03: " <<std::endl<<zk.topRows(3) <<std::endl;
	//std::cout << "zhat: " <<std::endl<<hnonlin2antenna(xBar,RR,Ls2p,Lcg2p).topRows(3)<<std::endl;
	Eigen::Matrix<double,6,15> Hk = getHkMatrixTwoAntennaTrueState(Ls2p,RR,Lcg2p);
	Eigen::Matrix<double,6,1> zMinusZbar = zk-hnonlin2antenna(xBar,RR,Ls2p,Lcg2p);
	Eigen::Matrix<double,6,6> Sk = Rbar+Hk*Pbar*Hk.transpose();
	Eigen::Matrix<double,15,6> Wk = Pbar*Hk.transpose()*Sk.inverse();
	xkp1 = xBar + Wk*zMinusZbar;
	//std::cout << "xkp1" << std::endl << xkp1 <<std::endl;
	Pkp1 = (Eigen::MatrixXd::Identity(15,15)-Wk*Hk)*Pbar*(Eigen::MatrixXd::Identity(15,15)-Wk*Hk).transpose()
		+ Wk*Rbar*Wk.transpose();
}


//Dynamic nonlinear propagation for IMU data. NOTE: THIS WANTS THE ORIGINAL fB, wB WITHOUT GRAVITY OR ROTMAT CORRECTION
Eigen::Matrix<double,15,1> gpsImuNode::fdyn(const Eigen::Matrix<double,15,1> x0, const double dt,
	const Eigen::Vector3d fB0, const Eigen::Vector3d wB0, const Eigen::Matrix3d RR, const Eigen::Vector3d lAB)
{
	Eigen::Matrix<double,15,1> x1 = x0;
	//split x0 into component vectors, assuming [p,v,gamma,ba,bg]
	const Eigen::Vector3d x = x0.topRows(3);
	const Eigen::Vector3d v = x0.middleRows(3,3);
	const Eigen::Vector3d gamma = x0.middleRows(6,3);
	const Eigen::Vector3d ba = x0.middleRows(9,3);
	const Eigen::Vector3d bg = x0.bottomRows(3);

	//Approximate propagation
	Eigen::Vector3d xkp1 = x + dt*v;
	Eigen::Vector3d omegaB = wB0-bg;
	Eigen::Vector3d wB_x_wB_x_lAB = omegaB.cross(omegaB.cross(lAB));
	Eigen::Vector3d a = RR.transpose()*(fB0 - wB_x_wB_x_lAB - ba) - Eigen::Vector3d(0,0,9.81);
	Eigen::Vector3d vkp1 = v + dt*a;
	Eigen::Vector3d gammakp1 = gamma + dt*omegaB;

	//Outputs--it is assumed that biases do not vary (time constant sufficiently long such that bkp1=bk)
	x1.topRows(3)=xkp1;
	x1.middleRows(3,3)=vkp1;
	x1.middleRows(6,3)=gammakp1;

	return x1;
}


//True state nonlinear measurement equation
Eigen::Matrix<double,6,1> gpsImuNode::hnonlin2antenna(const Eigen::Matrix<double,15,1> x0,
	const Eigen::Matrix3d RR, const Eigen::Vector3d ls2p, const Eigen::Vector3d lcg2p)
{
	Eigen::Matrix<double,6,1> zhat;

	Eigen::Matrix3d R2 = updateRBIfromGamma(RR,x0.middleRows(6,3));
	Eigen::Vector3d rCB = R2.transpose()*unit3(ls2p);
	zhat.topRows(3)=x0.topRows(3)+RR.transpose()*lcg2p;
	zhat.bottomRows(3)=rCB;
	return zhat;
}


//rotate by hatmat(gammavec) to new RBI
Eigen::Matrix3d gpsImuNode::updateRBIfromGamma(const Eigen::Matrix3d R0, const Eigen::Vector3d gamma)
{
    return (Eigen::Matrix3d::Identity()+hatmat(gamma))*R0;
}


//Generates F matrix from imu data
Eigen::Matrix<double,15,15> gpsImuNode::getFmatrixCF(const double dt, const Eigen::Vector3d fB,
	const Eigen::Vector3d omegaB, const Eigen::Matrix3d RR)
{
	//A is continuous time, Fk is discrete time
	Eigen::Matrix<double,15,15> Ak = Eigen::MatrixXd::Zero(15,15);
	Ak.block(0,3,3,3)=Eigen::Matrix3d::Identity();
	//Take off gravity here
	//Ak.block(3,6,3,3)=-RR.transpose()*hatmat(fB+RR.transpose()*Eigen::Vector3d(0,0,-9.81));
	//Measurement callback has already taken off gravity
	//Ak.block(3,6,3,3) = RR.transpose()*hatmat(fB-RR*Eigen::Vector3d(0,0,9.81));
	Ak.block(3,6,3,3) = RR.transpose()*hatmat(fB);
	Ak.block(3,9,3,3) = -RR.transpose();
	Ak.block(6,6,3,3) = RR*hatmat(omegaB);
	Ak.block(6,12,3,3) = -Eigen::MatrixXd::Identity(3,3);
	return Eigen::Matrix<double,15,15>::Identity() + dt*Ak;
}


//Grabs noise matrix
Eigen::Matrix<double,15,6> gpsImuNode::getGammakmatrixCF(const double dt, const Eigen::Matrix3d RR)
{
	Eigen::Matrix<double,15,6> Gammak = Eigen::Matrix<double,15,6>::Zero();
	Eigen::Matrix3d eye3 = Eigen::Matrix3d::Identity();
	Gammak.block(6,3,3,3)=1.0*dt*eye3;
	Gammak.block(0,0,3,3)=0.5*dt*dt*RR.transpose();
	Gammak.block(3,0,3,3)=1.0*dt*RR.transpose();
	Gammak.block(12,3,3,3)=eye3;
	Gammak.block(9,0,3,3)=eye3;
	return Gammak;
}


//Gets measurement equation for one antenna
Eigen::Matrix<double,3,15> gpsImuNode::getHkmatrixOneAntennaCF(const Eigen::Vector3d Lab, const Eigen::Matrix3d RR)
{
	Eigen::Matrix<double,3,15> Hk = Eigen::Matrix<double,3,15>::Zero();
	Hk.block(0,0,3,3)=Eigen::Matrix3d::Identity();
	Hk.block(0,6,3,3)=-RR*hatmat(Lab);
	return Hk;
}


//Gets measurement equation for two antennas
Eigen::Matrix<double,6,15> gpsImuNode::getHkmatrixTwoAntennaCF(const Eigen::Vector3d Limu,
	const Eigen::Vector3d Ls2p, const Eigen::Matrix3d RR)
{
	Eigen::Matrix<double,6,15> Hk = Eigen::Matrix<double,6,15>::Zero();
	Hk.block(0,0,3,3)=Eigen::Matrix3d::Identity();
	Hk.block(3,6,3,3)=RR.transpose()*hatmat(Ls2p);
	return Hk;
}


Eigen::Matrix<double,6,15> gpsImuNode::getHkMatrixTwoAntennaTrueState(const Eigen::Vector3d Ls2p,
	const Eigen::Matrix3d RR, const::Eigen::Vector3d Lcg2p)
{
	Eigen::Matrix<double,6,15> Hk = Eigen::Matrix<double,6,15>::Zero();
	Hk.block(0,0,3,3)=Eigen::Matrix3d::Identity();
	Hk.block(0,6,3,3)=RR.transpose()*hatmat(Lcg2p);
	Hk.block(3,6,3,3)=RR.transpose()*hatmat(unit3(Ls2p));
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


//Wabha solver.  Expects vI, vB as nx3 matrices with n sample vectors
Eigen::Matrix3d gpsImuNode::rotMatFromWahba(const Eigen::VectorXd weights,
	const::Eigen::MatrixXd vI, const::Eigen::MatrixXd vB)
{
	int n=weights.size();
	Eigen::Matrix3d B=Eigen::Matrix3d::Zero();
	for(int ij=0; ij<n; ij++)
	{
		B=B+weights(ij)*(vB.row(ij)).transpose()*vI.row(ij);
	}
	Eigen::Matrix3d U,V;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(B, Eigen::ComputeFullV | Eigen::ComputeFullU);
	U=svd.matrixU();
	V=svd.matrixV();
	Eigen::DiagonalMatrix<double, 3> M(1, 1, U.determinant()*V.determinant());
	//M.asDiagonal(Eigen::Vector3d(1,1,U.determinant()*V.determinant()));
	return U*M*(V.transpose());
}


Eigen::Vector3d gpsImuNode::unit3(const Eigen::Vector3d v1)
{
	return v1/v1.norm();
}


Eigen::Matrix3d gpsImuNode::orthonormalize(Eigen::Matrix3d)
{
	Eigen::Matrix3d outmat;
	
	
	return outmat;
}


//Adapted with minor changes from
//http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
Eigen::Quaterniond gpsImuNode::rotmat2quat(const Eigen::Matrix3d RR)
{
	double trace = RR.trace();
	Eigen::Matrix<double,4,1> q;
	if (trace > 0) {// RR_EPSILON = 0
		double s = 0.5 / sqrt(trace + 1.0);
		q << 0.25 / s,
			(RR(2,1) - RR(1,2))*s,
			(RR(0,2) - RR(2,0))*s,
			(RR(1,0) - RR(0,1))*s;
	}
	else {
		if (RR(0,0) > RR(1,1) && RR(0,0) > RR(2,2)) {
			double s = 2.0*sqrt(1.0 + RR(0,0) - RR(1,1) - RR(2,2));
			q << (RR(2,1) - RR(1,2))/s,
				 0.25*s,
				 (RR(0,1) + RR(1,0))/s,
				 (RR(0,2) + RR(2,0))/s;
		}
		else if (RR(1,1) > RR(2,2)) {
			double s = 2.0*sqrt(1.0 + RR(1,1) - RR(0,0) - RR(2,2));
			q << (RR(0,2) - RR(2,0))/s,
			     (RR(0,1) + RR(1,0))/s,
			     0.25 * s,
			     (RR(1,2) + RR(2,1))/s;
		}
		else {
			double s = 2.0*sqrt(1.0 + RR(2,2) - RR(0,0) - RR(1,1));
			q << (RR(1,0) - RR(0,1))/s,
			     (RR(0,2) + RR(2,0))/s,
			     (RR(1,2) + RR(2,1))/s,
			     0.25 * s;
		}
	}

	Eigen::Quaterniond quat;
	q=q/q.norm();
	quat.x()=q(1);
	quat.y()=q(2);
	quat.z()=q(3);
	quat.w()=q(0);
    return quat;	
}


Eigen::Matrix<double,15,15> gpsImuNode::getNumderivF(const double dv, const double dt,
	const Eigen::Matrix<double,15,1> x0,const Eigen::Vector3d fB0, const Eigen::Vector3d wB0,
	const Eigen::Matrix3d RR, const Eigen::Vector3d lAB)
{
	Eigen::Matrix<double,15,1> dForward, dBackward;
	Eigen::Matrix<double,15,15> FF;
	Eigen::Matrix<double,15,1> xPermuted;
	//Eigen::Matrix<double,15,1> xkp1=fdyn(x0,dt,fB0,wB0,RR,lAB);
	for(int ij=0; ij<15; ij++)
	{
		xPermuted=x0;
		xPermuted(ij)=xPermuted(ij)+dv;
		dForward=fdyn(xPermuted,dt,fB0,wB0,RR,lAB);
		xPermuted(ij)=xPermuted(ij)-2*dv;
		dBackward=fdyn(xPermuted,dt,fB0,wB0,RR,lAB);
//		FF.col(ij)=(dForward-xkp1)/dv;
		FF.col(ij)=(dForward-dBackward)/double(dv*2);
	}
	std::cout << "F_est" << std::endl << FF <<std::endl;
	return FF;
}

}

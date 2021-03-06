#include "filterImu.hpp"
#include <Eigen/SVD>


namespace gpsimu_odom
{
void gpsImu::runUKF(const imuMeas &imu, const gpsMeas &gps)
{
	double timu, tgps, dt0;
	Eigen::Vector3d imuAccelMeas, imuAttRateMeas, internal_rI, internal_rC;
	imu.getMeas(timu, imuAccelMeas, imuAttRateMeas); //Store current imu time in dt
	gps.getMeas(tgps, internal_rI, internal_rC);
	dt0=tgps-timu; //subtract off t0 to get actual dt
	
	//KF
	Eigen::Matrix<double,15,1> x0;
	Eigen::Matrix<double,15,15> P0;
	spkfPropagate15(xState_, Pimu_, (Qk12dividedByDt_*dt0), dt0, imuAccelMeas, imuAttRateMeas, RBI_, Limu_, x0, P0);
	updateRBIfromGamma(RBI_,x0);
	spkfMeasure6(x0, P0, Rk_, internal_rI, internal_rC, RBI_, Lcg2p_, Ls2p_, xState_, Pimu_);
	updateRBIfromGamma(RBI_,xState_);


	//Test bias saturation to avoid OOM errors
	//saturateBiases(maxBa_,maxBg_);

	RBI_=orthonormalize(RBI_);

	return;
}


void gpsImu::runUKFpropagateOnly(const double tPrev, const imuMeas &imu)
{
	double timu, dt0;
	Eigen::Vector3d imuAccelMeas, imuAttRateMeas;
	imu.getMeas(timu,imuAccelMeas,imuAttRateMeas); 
	dt0=timu-tPrev; //subtract off t0 to get actual dt
	//Qk12_ ~ (QK/dt_est)*dt
	Eigen::Matrix<double,15,1> x0;
	Eigen::Matrix<double,15,15> P0;
	spkfPropagate15(xState_,Pimu_,(Qk12dividedByDt_*dt0),dt0,imuAccelMeas,imuAttRateMeas,RBI_,Limu_, x0, P0);
	updateRBIfromGamma(RBI_,x0);

	xState_=x0;
	Pimu_=P0;

	//Test bias saturation to avoid OOM errors
	//saturateBiases(maxBa_,maxBg_);

	RBI_=orthonormalize(RBI_);

	return;
}


//Dynamic nonlinear propagation for IMU data. NOTE: This handles noise and noise rate for gyro/accel
Eigen::Matrix<double,15,1> gpsImu::fdynSPKF(const Eigen::Matrix<double,15,1> &x0, const double dt,
	const Eigen::Vector3d &fB0, const Eigen::Matrix<double,12,1> &vk, const Eigen::Vector3d &wB0,
	const Eigen::Matrix3d &RR, const Eigen::Vector3d &lAB)
{
	Eigen::Matrix<double,15,1> x1 = x0;
	//split x0 into component vectors, assuming [x,v,gamma,ba,bg]
	const Eigen::Vector3d x = x0.topRows(3);
	const Eigen::Vector3d v = x0.middleRows(3,3);
	const Eigen::Vector3d gamma = x0.middleRows(6,3);
	const Eigen::Vector3d ba = x0.middleRows(9,3);
	const Eigen::Vector3d bg = x0.bottomRows(3);
	const Eigen::Vector3d vgk = vk.topRows(3);
	const Eigen::Vector3d vgk2 = vk.middleRows(3,3);
	const Eigen::Vector3d vak = vk.middleRows(6,3);
	const Eigen::Vector3d vak2 = vk.bottomRows(3);
	const Eigen::Matrix3d RR2 = updateRBIfromGamma(RR,x0.middleRows(6,3));

	//Approximate propagation
	Eigen::Vector3d xkp1 = x + dt*v;
	Eigen::Vector3d omegaB = wB0 - bg - vgk;
	Eigen::Vector3d wB_x_wB_x_lAB = omegaB.cross(omegaB.cross(lAB));
	Eigen::Vector3d a = RR2.transpose()*(fB0 - wB_x_wB_x_lAB - ba - vak) - Eigen::Vector3d(0,0,9.8);
	Eigen::Vector3d vkp1 = v + dt*a;
	Eigen::Vector3d gammakp1 = gamma + dt*omegaB;

	//Outputs--it is assumed that biases do not vary (time constant sufficiently large such that bkp1=bk)
	x1.topRows(3) = xkp1;
	x1.middleRows(3,3) = vkp1;
	x1.middleRows(6,3) = gammakp1;
	x1.middleRows(9,3) = exp(-dt/tauA_)*x1.middleRows(9,3) + vak2;
	x1.bottomRows(3) = exp(-dt/tauG_)*x1.bottomRows(3) + vgk2;

	return x1;
}

//True state nonlinear measurement equation. lcg2p assumed a unit3 already
Eigen::Matrix<double,6,1> gpsImu::hnonlinSPKF(const Eigen::Matrix<double,15,1> &x0,
	const Eigen::Matrix3d &RR, const Eigen::Vector3d &ls2p, const Eigen::Vector3d &lcg2p,
	const Eigen::Matrix<double,6,1> &vk)
{
	Eigen::Matrix<double,6,1> zhat;

	Eigen::Matrix3d R2 = updateRBIfromGamma(RR,x0.middleRows(6,3));
	Eigen::Vector3d rCB = R2.transpose()*unit3(ls2p)+vk.bottomRows(3);
	zhat.topRows(3)=x0.topRows(3)+R2.transpose()*lcg2p+vk.topRows(3);
	zhat.bottomRows(3)=rCB;
	return zhat;
}


//Hardcoding matrix sizes instead of doing dynamic resizing to preserve speed
void gpsImu::spkfPropagate15(const Eigen::Matrix<double,15,1> &x0, const Eigen::Matrix<double,15,15> &P0,
	const Eigen::Matrix<double,12,12> &Q, const double dt, const Eigen::Vector3d &fB0, const Eigen::Vector3d &wB0,
	const Eigen::Matrix3d &RR, const Eigen::Vector3d &lAB, Eigen::Matrix<double,15,1> &xkp1, Eigen::Matrix<double,15,15> &Pkp1)
{
	static const double epsilon(1.0e-8);
	static const double alpha(1.0e-3);
	static const double beta(2.0);
	static const double kappa(0.0);
	static const int nn(27);
	static const double lambda = (alpha*alpha)*(kappa+double(nn))-double(nn);
	static const double w_mean_center(lambda/(double(nn)+lambda));
	static const double w_mean_reg(1.0/2.0/(double(nn)+lambda));
	static const double w_cov_center(w_mean_center+1.0-alpha*alpha+beta);
	static const double w_cov_reg(w_mean_reg);
	static const double cp(sqrt(double(nn)+lambda));
	Eigen::Matrix<double,27,27> Paug=Eigen::Matrix<double,27,27>::Zero();
	Eigen::Matrix<double,15,1> xBar, storeDum, xSPoint;
	Eigen::Matrix<double,15,55> xStore;
	Eigen::Matrix<double,27,27> cholP; // compute the Cholesky decomposition of A
	Paug.topLeftCorner(15,15)=P0;
	Paug.bottomRightCorner(12,12)=Q;

	cholP = (Paug.llt().matrixL());

	xBar = fdynSPKF(x0, dt, fB0, Eigen::Matrix<double,12,1>::Zero(), wB0, RR, lAB);
	xStore.col(0) = xBar;
	Eigen::Matrix<double,27,1> xAug, x_sp;
	xAug.topRows(15) = x0;
	xAug.bottomRows(12)=Eigen::Matrix<double,12,1>::Zero();
	xBar = w_mean_center*xBar;
	
	int colno;
	double spSign; //This can be an int (will only have values on +-1), I'm just being careful to avoid implicit conversions.
	//Propagate through sigma points
	spSign=1.0;
	for(int ij=0; ij<2*nn; ij++)
	{
		//std::cout << "XBAR:" << std::endl << xBar <<std::endl;
		colno = ij%nn;
		if(ij>=nn)
		{
			spSign=-1.0;
		}
		x_sp = xAug + cp*spSign*cholP.col(colno);
		xSPoint=x_sp.topRows(15);
		storeDum = fdynSPKF(xSPoint, dt, fB0, x_sp.bottomRows(12), wB0, RR, lAB);
		xStore.col(ij+1) = storeDum;
		xBar = xBar + w_mean_reg*storeDum;
	}
	//std::cout << "xbar:" <<std::endl<<xBar<<std::endl;

	//Recombine for covariance
	Eigen::Matrix<double,15,15> PbaRk_p1;
	PbaRk_p1 = w_cov_center*(xStore.col(0)-xBar)*((xStore.col(0)-xBar).transpose());
	for(int ij=0; ij<2*nn; ij++)
	{
		PbaRk_p1 = PbaRk_p1 + w_cov_reg*(xStore.col(ij+1)-xBar)*((xStore.col(ij+1)-xBar).transpose());
	}
	//std::cout << "Pmax: " << PbaRk_p1.maxCoeff() << std::endl;

	//outputs
	xkp1 = xBar;
	Pkp1 = PbaRk_p1;
	//std::cout << "P:" <<std::endl << Pkp1 <<std::endl;
	return;
}


//Hardcoding matrix sizes instead of doing dynamic resizing to preserve speed
//lcg2p, ls2p are the real values of the antenna locations in the body.  ls2p will be unit3'd inside of this function.
void gpsImu::spkfMeasure6(const Eigen::Matrix<double,15,1> &x0, const Eigen::Matrix<double,15,15> &P0,
	const Eigen::Matrix<double,6,6> &R, const Eigen::Vector3d &rI_measurement, const Eigen::Vector3d &rCu_measurement,
	const Eigen::Matrix3d &RR, const Eigen::Vector3d &lcg2p, const Eigen::Vector3d &ls2p,
	Eigen::Matrix<double,15,1> &xkp1, Eigen::Matrix<double,15,15> &Pkp1)
{
	static const double epsilon(1.0e-8);
	static const double alpha(1.0e-3);
	static const double beta(2.0);
	static const double kappa(0.0);
	static const int nn(21);
	static const double lambda = (alpha*alpha)*(kappa+double(nn))-double(nn);
	static const double w_mean_center(lambda/(double(nn)+lambda));
	static const double w_mean_reg(1.0/2.0/(double(nn)+lambda));
	static const double w_cov_center(w_mean_center+1.0-alpha*alpha+beta);
	static const double w_cov_reg(w_mean_reg);
	static const double cp(sqrt(nn+lambda));
	Eigen::Matrix<double,21,21> Paug=Eigen::Matrix<double,21,21>::Zero();
	Eigen::Matrix<double,6,1> zBar, storeDum;
	Eigen::Matrix<double,6,43> zStore;
	Eigen::Matrix<double,15,43> xStore;
	Eigen::Matrix<double,15,1> xSPoint;
	Eigen::Matrix<double,21,21> cholP; // compute the Cholesky decomposition of A
	Eigen::Vector3d ls2p_unit3 = unit3(ls2p);

	//Center point and propagation
	Paug.topLeftCorner(15,15)=P0; Paug.bottomRightCorner(6,6)=R;
	cholP = (Paug.llt().matrixL());
	zBar = hnonlinSPKF(x0, RR, ls2p, lcg2p, Eigen::Matrix<double,6,1>::Zero());
	zStore.col(0) = zBar;
	xStore.col(0) = x0;
	Eigen::Matrix<double,21,1> xAug, x_sp;
	xAug.topRows(15) = x0;
	xAug.bottomRows(6)=Eigen::Matrix<double,6,1>::Zero();
	zBar = w_mean_center*zBar;
	int colno;
	double spSign; //This can be an int (will only have values on +-1), I'm just being careful to avoid implicit conversions.

	//Propagate regression points
	spSign=1.0;
	for(int ij=0; ij<2*nn; ij++)
	{
		colno = ij%nn;
		if(ij>=nn)
		{
			spSign=-1.0;
		}
		x_sp = xAug + cp*spSign*cholP.col(colno);
		xSPoint = x_sp.topRows(15);
		storeDum = hnonlinSPKF(xSPoint, RR, ls2p, lcg2p, x_sp.bottomRows(6));
		zStore.col(ij+1) = storeDum;
		xStore.col(ij+1) = x_sp.topRows(15);
		zBar = zBar + w_mean_reg*storeDum;
	}

	//Recombine for covariance
	Eigen::Matrix<double,15,6> Pxz;
	Eigen::Matrix<double,6,6> Pzz;
	Eigen::Matrix<double,6,1> dz;
	Pxz = w_cov_center*(xStore.col(0)-x0)*((zStore.col(0)-zBar).transpose());
	Pzz = w_cov_center*(zStore.col(0)-zBar)*((zStore.col(0)-zBar).transpose());
	for(int ij=0; ij<2*nn; ij++)
	{
		dz = zStore.col(ij+1)-zBar; //For efficiency
		Pxz = Pxz + w_cov_reg*(xStore.col(ij+1)-x0)*(dz.transpose());
		Pzz = Pzz + w_cov_reg*dz*(dz.transpose());
	}

	//LMMSE
	Eigen::Matrix<double,6,1> z_measurement;
	const Eigen::Matrix<double,6,6> PzzInv = Pzz.inverse(); //Efficiency
	z_measurement.topRows(3)=rI_measurement;
	z_measurement.bottomRows(3)=unit3(rCu_measurement);
	xkp1 = x0 + Pxz*PzzInv*(z_measurement-zBar);
	Pkp1 = P0 - Pxz*PzzInv*(Pxz.transpose());
	return;
}


//Rotation matrix
Eigen::Matrix3d gpsImu::euler2dcm312(const Eigen::Vector3d &ee)
{
  	const double cPhi = cos(ee(0));
  	const double sPhi = sin(ee(0));
  	const double cThe = cos(ee(1));
  	const double sThe = sin(ee(1));
  	const double cPsi = cos(ee(2));
  	const double sPsi = sin(ee(2));
  	Eigen::Matrix3d R2;
  	R2 << cPsi*cThe - sPhi*sPsi*sThe, cThe*sPsi + cPsi*sPhi*sThe, -cPhi*sThe,
        -cPhi*sPsi,                                  cPhi*cPsi,       sPhi,
        cPsi*sThe + cThe*sPhi*sPsi, sPsi*sThe - cPsi*cThe*sPhi,  cPhi*cThe;
  	return R2;
}


//rotate by hatmat(gammavec) to new RBI_
Eigen::Matrix3d gpsImu::updateRBIfromGamma(const Eigen::Matrix3d &R0, const Eigen::Vector3d &gamma)
{
    //return (Eigen::Matrix3d::Identity()+hatmat(gamma))*R0;
    //return (rotMatFromEuler(gamma))*R0;
    return (euler2dcm312(gamma))*R0;
}


//rotate by hatmat(gammavec) to new RBI_
void gpsImu::updateRBIfromGamma(Eigen::Matrix3d &R0, Eigen::Matrix<double,15,1> &x0)
{
	R0=euler2dcm312(x0.middleRows(6,3))*R0;
	x0.middleRows(6,3)=Eigen::Vector3d::Zero();
}


//Cross product equivalent.  Named this way for consistency with nn_imu_dat
Eigen::Matrix3d gpsImu::hatmat(const Eigen::Vector3d &v1)
{
	Eigen::Matrix3d f = Eigen::MatrixXd::Zero(3,3);
	f(0,1)=-v1(2); f(0,2)=v1(1);
	f(1,0)=v1(2); f(1,2)=-v1(0);
	f(2,0)=-v1(1); f(2,1)=v1(0);
	return f;
}


//Wabha solver.  Expects vI, vB as nx3 matrices with n sample vectors
Eigen::Matrix3d gpsImu::rotMatFromWahba(const Eigen::VectorXd &weights,
	const::Eigen::MatrixXd &vI, const::Eigen::MatrixXd &vB)
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


Eigen::Vector3d gpsImu::unit3(const Eigen::Vector3d &v1)
{
	return v1/v1.norm();
}


double gpsImu::symmetricSaturationDouble(const double inval, const double maxval)
{
	//Handling this with an error to avoid needing an abs() each call
	if(maxval<0)
	{ std::cout <<"ERROR: Saturation bound is negative" << std::endl;}
	if(inval > maxval)
	{
		return maxval;
	}else if(inval < -1.0*maxval)
	{
		return -1.0*maxval;
	}else
	{
		return inval;
	}
}


void gpsImu::saturateBiases(const double baMax, const double bgMax)
{
	//Saturation
	for(int ij=0; ij<3; ij++)
	{
		xState_(ij+9) = symmetricSaturationDouble(xState_(ij+9),baMax);
	}
	for(int ij=0; ij<3; ij++)
	{
		xState_(ij+12) = symmetricSaturationDouble(xState_(ij+12),bgMax);
	}
}



Eigen::Matrix3d gpsImu::orthonormalize(const Eigen::Matrix3d &inmat)
{
	Eigen::Matrix3d outmat;
	Eigen::Matrix3d U,V;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(inmat, Eigen::ComputeFullV | Eigen::ComputeFullU);
	U=svd.matrixU();
	V=svd.matrixV();
	outmat = U*(V.transpose());
	return outmat;
}


//Adapted with minor changes from
//http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
Eigen::Quaterniond gpsImu::rotmat2quat(const Eigen::Matrix3d &RR)
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


Eigen::Matrix3d gpsImu::rotMatFromQuat(const Eigen::Quaterniond &qq)
{
	const double xx=qq.x();
	const double yy=qq.y();
	const double zz=qq.z();
	const double ww=qq.w();
	Eigen::Matrix3d RR;
  /*RR << 1-2*yy*yy-2*zz*zz, 2*xx*yy+2*ww*zz, 2*xx*zz-2*ww*yy,
        2*xx*yy-2*ww*zz, 1-2*xx*xx-2*zz*zz, 2*yy*zz+2*ww*xx,
        2*xx*zz+2*ww*yy, 2*yy*zz-2*ww*xx, 1-2*xx*xx-2*yy*yy; // the transposed derp*/
	RR << 1-2*yy*yy-2*zz*zz, 2*xx*yy-2*ww*zz, 2*xx*zz+2*ww*yy,
    	2*xx*yy+2*ww*zz, 1-2*xx*xx-2*zz*zz, 2*yy*zz-2*ww*xx,
    	2*xx*zz-2*ww*yy, 2*yy*zz+2*ww*xx, 1-2*xx*xx-2*yy*yy;
	return RR;
}


Eigen::Matrix3d gpsImu::rotMatFromEuler(const Eigen::Vector3d &ee)
{
  const double phi=ee(0);
  const double theta=ee(1);
  const double psi=ee(2);
  Eigen::Matrix3d RR;
  RR<<cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta),
      sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), sin(theta)*sin(phi)*sin(psi)+cos(phi)*cos(psi), sin(phi)*cos(theta),
      cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), cos(phi)*cos(theta);
  return RR;
}

}

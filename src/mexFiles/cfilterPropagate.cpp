#include <Eigen/Eigen>
#include "mex.h"


//Function prototypes
void unpackMxArray(const double * y, const int mRows, const int nCols, Eigen::MatrixXd &outmat);
//You can add more prototypes if your code has multiple functions.  PUT YOUR CODE HERE
void spkfPropagate15(const Eigen::Matrix<double,15,1> x0, const Eigen::Matrix<double,15,15> P0,
	const Eigen::Matrix<double,12,12> Q, const double dt, const Eigen::Vector3d fB0, const Eigen::Vector3d wB0,
	const Eigen::Matrix3d RR, const Eigen::Vector3d lAB, Eigen::Matrix<double,15,15> &Pkp1, Eigen::Matrix<double,15,1> &xkp1);

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	//input order: xk, Pk, Qk, dt, RBI, fB, wB, L_imu
	//INPUTS
 	//it is assumed that the user will enter at most 16 inputs.  This can be changed.
 	int numberOfInputs = 8;
 	Eigen::MatrixXd datastore[numberOfInputs];
 	for(int ij=0; ij<numberOfInputs; ij++)
 	{
 		double *y = mxGetPr(prhs[ij]);
  		int mrows = mxGetM(prhs[ij]);
 		int ncols = mxGetN(prhs[ij]); 	
 		unpackMxArray(mxGetPr(prhs[ij]),mrows,ncols,datastore[ij]);
 	}

 	//Move elements from datastore to Eigen matrices
 	Eigen::MatrixXd xk,Pk,Qk,dt,RBI,fB,wB,L_imu;
 	xk = datastore[0];
 	Pk = datastore[1];
 	Qk = datastore[2];
 	dt = datastore[3];
 	RBI = datastore[4];
 	fB = datastore[5];
 	wB = datastore[6];
 	L_imu = datastore[7];

 	//PROCESSING
 	Eigen::MatrixXd testOut;
 	Eigen::Matrix<double,15,1> xkp1;
 	Eigen::Matrix<double,15,15> Pkp1;
 	spkfPropagate15(xk, Pk, Qk, dt(0), fB, wB, RBI, L_imu, Pkp1, xkp1);


 	//OUTPUTS
 	int numoutputs=2;
 	Eigen::MatrixXd outputs[numoutputs];
 	//Determine which variables to output
 	outputs[0] = xkp1;
 	outputs[1] = Pkp1;

 	//Package elements from outputs[] into the desired
 	mwSize dims[2];
 	int row_max;
 	for(int ij=0; ij < numoutputs; ij++)
    {
		dims[0] = (unsigned int)outputs[ij].rows();
		dims[1] = (unsigned int)outputs[ij].cols();
		plhs[ij] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
		double *outMatrix_ptr = mxGetPr(plhs[ij]);
    	for (int col=0; col < outputs[ij].cols(); col++)
   	 	{
    		row_max = outputs[ij].rows();
	        for (int row=0; row < row_max; row++)
        	{
            	outMatrix_ptr[row + col*row_max] = outputs[ij](row,col);
        	}
    	}
	}

	//Don't forget the return; if you do, Matlab will crash.
    return;
}


void unpackMxArray(const double * y, const int mRows, const int nCols, Eigen::MatrixXd &outmat)
{
	Eigen::MatrixXd output;
	output.resize(mRows,nCols);
	for(int iRow=0; iRow<mRows; iRow++)
	{
		for(int iCol=0; iCol<nCols; iCol++)
		{
			output(iRow,iCol) = y[iRow+iCol*mRows];
		}
	}
	outmat = output;
}

//Cross product equivalent.  Named this way for consistency with nn_imu_dat
Eigen::Matrix3d hatmat(const Eigen::Vector3d v1)
{
	Eigen::Matrix3d f = Eigen::MatrixXd::Zero(3,3);
	f(0,1)=-v1(2); f(0,2)=v1(1);
	f(1,0)=v1(2); f(1,2)=-v1(0);
	f(2,0)=-v1(1); f(2,1)=v1(0);
	return f;
}



Eigen::Matrix3d rotMatFromEuler(Eigen::Vector3d ee)
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


//Rotation matrix
Eigen::Matrix3d euler2dcm312(const Eigen::Vector3d ee)
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


//rotate by hatmat(gammavec) to new RBI
Eigen::Matrix3d updateRBIfromGamma(const Eigen::Matrix3d R0, const Eigen::Vector3d gamma)
{
    //return (Eigen::Matrix3d::Identity()+hatmat(gamma))*R0;
    //return (rotMatFromEuler(gamma))*R0;
    return (euler2dcm(gamma))*R0;
}

Eigen::Vector3d unit3(const Eigen::Vector3d v1)
{
	return v1/v1.norm();
}

//Dynamic nonlinear propagation for IMU data. NOTE: This handles noise and noise rate for gyro/accel
Eigen::Matrix<double,15,1> fdynSPKF(const Eigen::Matrix<double,15,1> x0, const double dt,
	const Eigen::Vector3d fB0, const Eigen::Matrix<double,12,1> vk, const Eigen::Vector3d wB0,
	const Eigen::Matrix3d RR, const Eigen::Vector3d lAB)
{
	static const double tauA(100.0), tauG(100.0);
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
	Eigen::Vector3d a = RR2.transpose()*(fB0 - ba - wB_x_wB_x_lAB - vak) - Eigen::Vector3d(0,0,9.8);
	Eigen::Vector3d vkp1 = v + dt*a;
	Eigen::Vector3d gammakp1 = gamma + dt*omegaB;

	//Outputs--it is assumed that biases do not vary (time constant sufficiently large such that bkp1=bk)
	x1.topRows(3) = xkp1;
	x1.middleRows(3,3) = vkp1;
	x1.middleRows(6,3) = gammakp1;
	x1.middleRows(9,3) = exp(-dt/tauA)*x1.middleRows(9,3) + vak2;
	x1.bottomRows(3) = exp(-dt/tauG)*x1.bottomRows(3) + vgk2;

	return x1;
}

//True state nonlinear measurement equation. lcg2p assumed a unit3 already
Eigen::Matrix<double,6,1> hnonlinSPKF(const Eigen::Matrix<double,15,1> x0,
	const Eigen::Matrix3d RR, const Eigen::Vector3d ls2p, const Eigen::Vector3d lcg2p,
	const Eigen::Matrix<double,6,1> vk)
{
	Eigen::Matrix<double,6,1> zhat;

	Eigen::Matrix3d R2 = updateRBIfromGamma(RR,x0.middleRows(6,3));
	Eigen::Vector3d rCB = R2.transpose()*unit3(ls2p+vk.bottomRows(3));
	zhat.topRows(3)=x0.topRows(3)+R2.transpose()*lcg2p+vk.topRows(3);
	zhat.bottomRows(3)=rCB;
	return zhat;
}


//Hardcoding matrix sizes instead of doing dynamic resizing to preserve speed
void spkfPropagate15(const Eigen::Matrix<double,15,1> x0, const Eigen::Matrix<double,15,15> P0,
	const Eigen::Matrix<double,12,12> Q, const double dt, const Eigen::Vector3d fB0, const Eigen::Vector3d wB0,
	const Eigen::Matrix3d RR, const Eigen::Vector3d lAB, Eigen::Matrix<double,15,15> &Pkp1, Eigen::Matrix<double,15,1> &xkp1)
{
	//std::cout << "dt: " << dt <<std::endl;
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
	Eigen::Matrix<double,15,1> xBar, storeDum;
	Eigen::Matrix<double,15,55> xStore;
	Eigen::Matrix<double,27,27> cholP; // compute the Cholesky decomposition of A
	Paug.topLeftCorner(15,15)=P0;
	Paug.bottomRightCorner(12,12)=Q;
	//std::cout << "Eig(P):" << std::endl << Paug.eigenvalues() <<std::endl;

	cholP = (Paug.llt().matrixL());
	
	//std::cout << "chol(P_aug)" <<std::endl<<cholP<<std::endl; 
	//std::cout << "maxval: " << cholP.maxCoeff() << std::endl;
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

	Eigen::Matrix<double,27,54> sps = inTestMat;

	for(int ij=0; ij<2*nn; ij++)
	{
		//std::cout << "XBAR:" << std::endl << xBar <<std::endl;
		colno = ij%nn;
		if(ij>=nn)
		{
			colno = ij-nn;
			spSign=-1.0;
		}
		x_sp = xAug + cp*spSign*cholP.col(colno);
		storeDum = fdynSPKF(x_sp.topRows(15), dt, fB0, x_sp.bottomRows(12), wB0, RR, lAB);
		xStore.col(ij+1) = storeDum;
		xBar = xBar + w_mean_reg*storeDum;
	}
	//std::cout << "xbar:" <<std::endl<<xBar<<std::endl;

	//Recombine for covariance
	Eigen::Matrix<double,15,15> Pbarkp1;
	Pbarkp1 = w_cov_center*(xStore.col(0)-xBar)*((xStore.col(0)-xBar).transpose());
	for(int ij=0; ij<2*nn; ij++)
	{
		Pbarkp1 = Pbarkp1 + w_cov_reg*(xStore.col(ij+1)-xBar)*((xStore.col(ij+1)-xBar).transpose());
	}
	//std::cout << "Pmax: " << Pbarkp1.maxCoeff() << std::endl;

	outTestMat = xStore;

	//outputs
	xkp1 = xBar;
	Pkp1 = Pbarkp1;
	//std::cout << "P:" <<std::endl << Pkp1 <<std::endl;
	return;
}


//Hardcoding matrix sizes instead of doing dynamic resizing to preserve speed
//lcg2p, ls2p are the real values of the antenna locations in the body.  ls2p will be unit3'd inside of this function.
void spkfMeasure6(const Eigen::Matrix<double,15,1> x0, const Eigen::Matrix<double,15,15> P0,
	const Eigen::Matrix<double,6,6> R, const Eigen::Vector3d rI_measurement, const Eigen::Vector3d rCu_measurement,
	const Eigen::Matrix3d RR, const Eigen::Vector3d lcg2p, const Eigen::Vector3d ls2p,
	Eigen::Matrix<double,15,15> &Pkp1, Eigen::Matrix<double,15,1> &xkp1)
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
	Eigen::Matrix<double,21,21> cholP; // compute the Cholesky decomposition of A
	Eigen::Vector3d ls2p_unit3 = unit3(ls2p);

	//Center point and propagation
	Paug.topLeftCorner(15,15)=P0; Paug.bottomRightCorner(6,6)=R;
	cholP = (Paug.llt().matrixL()).transpose();
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
		storeDum = hnonlinSPKF(x_sp.topRows(15), RR, ls2p, lcg2p, x_sp.bottomRows(6));
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




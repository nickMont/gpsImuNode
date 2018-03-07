#include <iostream>
#include <eigen3/Eigen/Geometry>

class gnssimu
{
public:

    //constructor
	gnssimu();

	//working functions
	Eigen::Matrix3d euler_to_CTM(Eigen::Vector3d eul);
	Eigen::Vector3d CTM_to_euler(Eigen::Matrix3d CC);
	Eigen::Matrix3d hatmat(Eigen::Vector3d vec);
	Eigen::Vector3d ecef2lla(Eigen::Vector3d ecef);
	Eigen::Vector3d unit3(Eigen::Vector3d &x_in);

	//workhorses
	Eigen::VectorXd initialize_enu(Eigen::VectorXd init);
	void updateState(Eigen::VectorXd instate);
	Eigen::MatrixXd returnHcoupled();
	Eigen::MatrixXd returnFimu(double dt);
	Eigen::VectorXd returnMeasurement();
	void processIMU(Eigen::VectorXd measIMU);
	void processGPS(Eigen::VectorXd measGPS);

	//KF
	void kalmanMeasure(Eigen::VectorXd &xkbar, Eigen::MatrixXd &Pkbar, Eigen::VectorXd &zMeas,
		Eigen::MatrixXd &H, Eigen::MatrixXd &R, Eigen::VectorXd &xkp1, Eigen::MatrixXd &Pkp1);
	void kalmanPropagate(Eigen::VectorXd &xk, Eigen::MatrixXd &Pk, Eigen::MatrixXd &F,
		Eigen::MatrixXd &Q, Eigen::MatrixXd &Gamma, Eigen::VectorXd &xkbar, Eigen::MatrixXd &Pkbar);
	void gpsImuCallback(Eigen::VectorXd zMeas);
	void imuCallback(Eigen::VectorXd zMeas);
	void gpsOnlyCallback(Eigen::VectorXd zMeas);

	//getters and setters
	Eigen::Vector3d returnXecef() {return x_ecef;}
	Eigen::Vector3d returnVecef() {return v_ecef;}
	Eigen::VectorXd returnEVX() {return x_reported;}
	bool setInitStates(Eigen::VectorXd initstates, Eigen::Vector3d leverarm);
	bool setInitCovariance(Eigen::MatrixXd initCov, Eigen::MatrixXd initMeasCov);

private:
	Eigen::Vector3d lever_ab, euler, rIMU, rGPS, vIMU, vGPS, b_a, b_g, omega_hat, x_ecef, v_ecef;
	Eigen::VectorXd imustate, x_reported, evxState;  //x_reported is the output euler/vel/pose
	double tlastgps, tlastimu;
	Eigen::Matrix3d C_hat, Omega_mat;
	Eigen::MatrixXd Qimu, Rgps, Pimu, Fimu, Gamma;
	bool initImu, initGps;
};
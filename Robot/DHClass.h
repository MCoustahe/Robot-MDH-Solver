#pragma once

#include "iostream"
#include "fstream"
#include "Eigen/Dense"
#include "vector"

using namespace std;
using namespace Eigen;

#define	pi  3.14159265358979323846264338327950288419716939;	// The magic number
#define	deg_rad  57.295779513082;									// Degrees to radians factor

#define MAXSTRLEN  2048 /* 2K */
#define SKIP_LINE(f){                                                       \
	char buf[MAXSTRLEN];                                                        \
			while(!feof(f))                                                           \
	if(!fgets(buf, MAXSTRLEN-1, f) || buf[strlen(buf)-1]=='\n') break;      \
}


class DH
{
public:
	DH();
	DH(double a[6], double alpha[6],double d[6],double theta[6]);
	DH(bool key);
	
	void printlength_a() const;
	void printtwist_alpha() const;
	void printoffset_d() const;
	void printangle_theta() const;

	void printMatA() const;
	void printMatD() const;
	void printMatAlpha() const;
	void printMatTheta() const;
	void printMatJ() const;

	void CreateMatAWithError(vector<double>& stdA);
	void CreateMatDWithError(vector<double>& stdD);
	void CreateMatAlphaWithError(vector<double>& stdAlpha);
	void CreateMatThetaWithError(vector<double>& stdTheta);
	void CreateJointWithError(vector<vector<double>>& std);

	void CreateDiffJoint(vector<vector<double>>& std);
	int DH::UpdateDiffJoint(vector<double>& theta, vector<vector<double>>& std);

	vector<Matrix4d> CalEndPose(vector<vector<double>>& Theta);
	vector<Matrix4d> CalEndPoseWithError(vector<vector<double>>& Theta, vector<vector<double>>& std);
	
	vector<Vector4d> CalTargetPose(vector<Matrix4d>& EndPose, Vector4d b);
	vector<Vector3d> CalDiff_Tar2Laser(vector<Vector4d> TargetPose, vector<Vector4d> Laser, Matrix4d Trans_Tar2Laser);
	vector<Vector3d> CalDiff_Laser2Tar(vector<Vector4d> Laser, vector<Vector4d> TargetPose, Matrix4d Trans_Tar2Laser);
	
	MatrixXd TransDiff2Fp(vector<Vector3d>& TmpDiff);

	void SumDiff(vector<Vector3d> Diff);

	int MDH_Solver(double Sum_Delta, int MaxIter, double MinError, double MinDelta, 
					vector<vector<double>> Theta, vector<Vector4d> Laser, Vector4d b1, Matrix4d Trans_Tar2Laser);


private:
	double length_a[6];
	double twist_alpha[6];
	double offset_d[6];
	double angle_theta[6];

	bool ErrorModelFinishedBtn = false;

	std::vector<Eigen::Matrix4d> MatA, MatAlpha, MatD, MatTheta, MatJ;
	std::vector<Eigen::Matrix4d> MateA, MateAlpha, MateD, MateTheta, MateJ;
	std::vector<Eigen::Matrix4d> MatdA, MatdAlpha, MatdD, MatdTheta, MatdJ;

	vector<Vector3d> diff_all_T2L;
	vector<Vector3d> diff_all_L2T;
	double Sum_Diff;
	
	char* szDHpara = "../data/DHpara.txt";

	int ReadDHpara(FILE* fp);
	void CreateInitJoint();

	int UpdateJoint(std::vector<double>& theta);
	Matrix4d CalEndPose();
	int UpdateJointWithError(vector<double>& theta, vector<vector<double>>& std);
	Matrix4d CalEndPoseWithError();

	Matrix4d dAlp(int i);
	Matrix4d dA(int i);
	Matrix4d dThe(int i);
	Matrix4d dD(int i);

};




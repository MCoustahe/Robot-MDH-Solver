#include "DHClass.h"

using namespace std;


DH::DH():length_a{0,0,0,0,0,0}, twist_alpha{ 0,0,0,0,0,0 }, offset_d{ 0,0,0,0,0,0 }, angle_theta{ 0,0,0,0,0,0 }
{}

DH::DH(double a[6], double alpha[6], double d[6], double theta[6]) 
		:length_a{  a[0],
					a[1],
					a[2],
					a[3],
					a[4],
					a[5] },
		 twist_alpha{   alpha[0],
						alpha[1],
						alpha[2],
						alpha[3],
						alpha[4],
						alpha[5] },
		 offset_d{  d[0],
					d[1],
					d[2],
					d[3],
					d[4],
					d[5] },
		 angle_theta{	theta[0],
						theta[1],
						theta[2],
						theta[3],
						theta[4],
						theta[5] }
{}

DH::DH(bool key)
{
	if (key) {
		FILE* fpDHpara = fopen(szDHpara, "r");
		ReadDHpara(fpDHpara);

		CreateInitJoint();
	}
}

//////////////////////////////////////////
// Print Function
void DH::printlength_a() const
{
	for (size_t i = 0; i < 6; i++) 
	{
		cout << "length_a" << i + 1 << " : " << length_a[i] << " ";
		if (i == 2)
			cout << endl;
	}
	cout << endl;
}

void DH::printtwist_alpha() const
{
	for (size_t i = 0; i < 6; i++)
	{
		cout << "twist_alpha" << i + 1 << " : " << twist_alpha[i] << " ";
		if (i == 2)
			cout << endl;
	}
	cout << endl;
}

void DH::printoffset_d() const
{
	for (size_t i = 0; i < 6; i++)
	{
		cout << "offset_d" << i + 1 << " : " << offset_d[i] << " ";
		if (i == 2)
			cout << endl;
	}
	cout << endl;
}

void DH::printangle_theta() const
{
	for (size_t i = 0; i < 6; i++)
	{
		cout << "angle_theta" << i + 1 << " : " << angle_theta[i] << " ";
		if (i == 2)
			cout << endl;
	}
	cout << endl;
}

void DH::printMatA() const {
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 4; j++){
			for (int k = 0; k < 4; k++) cout << MatA[i](j, k) << " ";
			cout << endl;
		}
	}
}

void DH::printMatD() const {
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) cout << MatD[i](j, k) << " ";
			cout << endl;
		}
	}
}

void DH::printMatAlpha() const {
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) cout << MatAlpha[i](j, k) << " ";
			cout << endl;
		}
	}
}

void DH::printMatTheta() const {
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) cout << MatTheta[i](j, k) << " ";
			cout << endl;
		}
	}
}

void DH::printMatJ() const {
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) cout << MatJ[i](j, k) << " ";
			cout << endl;
		}
	}
}

/////////////////////////////
// Create Error Joint Matrix
/////////////////////////////

void DH::CreateMatAWithError(vector<double>& stdA) {
	if (stdA.size() == 6) {
		for (int i = 0; i < stdA.size(); i++)
		{
			Matrix4d A = MatrixXd::Identity(4, 4);
			A(0, 3) = length_a[i] + stdA[i];
			MateA.push_back(A);
		}
	}
}

void DH::CreateMatDWithError(vector<double>& stdD) {
	if (stdD.size() == 6) {
		for (int i = 0; i < stdD.size(); i++)
		{
			Matrix4d B = MatrixXd::Identity(4, 4);
			B(2, 3) = offset_d[i] + stdD[i];
			MateD.push_back(B);
		}
	}
}

void DH::CreateMatAlphaWithError(vector<double>& stdAlpha) {
	if (stdAlpha.size() == 6) {
		for (int i = 0; i < stdAlpha.size(); i++)
		{
			double tmpalp = twist_alpha[i] + stdAlpha[i];
			Matrix4d C = MatrixXd::Identity(4, 4);
			if (abs(cos(tmpalp)) < 0.0000001)
			{
				C(1, 1) = 0;
				C(1, 2) = -sin(tmpalp);
				C(2, 1) = sin(tmpalp);
				C(2, 2) = 0;
			}
			else if (abs(sin(tmpalp)) < 0.0000001)
			{
				C(1, 1) = cos(tmpalp);
				C(1, 2) = 0;
				C(2, 1) = 0;
				C(2, 2) = cos(tmpalp);
			}
			else
			{
				C(1, 1) = cos(tmpalp);
				C(1, 2) = -sin(tmpalp);
				C(2, 1) = sin(tmpalp);
				C(2, 2) = cos(tmpalp);
			}
			MateAlpha.push_back(C);
		}
	}
}

void DH::CreateMatThetaWithError(vector<double>& stdTheta) {
	if (stdTheta.size() == 6) {
		for (int i = 0; i < stdTheta.size(); i++)
		{
			double tmpthe = angle_theta[i] + stdTheta[i];
			Matrix4d D = MatrixXd::Identity(4, 4);

			if (abs(cos(tmpthe)) < 0.0000001)
			{
				D(0, 0) = 0;
				D(0, 1) = -sin(tmpthe);
				D(1, 0) = sin(tmpthe);
				D(1, 1) = 0;
			}
			else if (abs(sin(tmpthe)) < 0.0000001)
			{
				D(0, 0) = cos(tmpthe);
				D(0, 1) = 0;
				D(1, 0) = 0;
				D(1, 1) = cos(tmpthe);
			}
			else
			{
				D(0, 0) = cos(tmpthe);
				D(0, 1) = -sin(tmpthe);
				D(1, 0) = sin(tmpthe);
				D(1, 1) = cos(tmpthe);
			}
			MateTheta.push_back(D);
		}
	}
}

void DH::CreateJointWithError(vector<vector<double>>& std) {

	CreateMatAWithError(std[0]);
	CreateMatDWithError(std[1]);
	CreateMatAlphaWithError(std[2]);
	CreateMatThetaWithError(std[3]);

	if (MateA.size() == 6 && MateD.size() == 6 && MateAlpha.size() == 6 && MateTheta.size() == 6)
		for (int i = 0; i < 6; i++)
		{
			Matrix4d E = MateAlpha[i] * MateA[i] * MateTheta[i] * MateD[i];
			MateJ.push_back(E);
		}
	if (MateJ.size() == 6)
		ErrorModelFinishedBtn = true;
}

/////////////////////////////
// Create Differential Joint Matrix
/////////////////////////////

void DH::CreateDiffJoint(vector<vector<double>>& std) {
	vector<double> stdalp, stdthe;
	if (std.size() == 4) {
		stdalp = std[2];
		stdthe = std[3];
	}

	if (std[2].size() == 6) {
		for (int i = 0; i < std[2].size(); i++)
		{
			double tmpalp = twist_alpha[i] + stdalp[i];
			double tmpthe = angle_theta[i] + stdthe[i];

			Matrix4d C = MatrixXd::Zero(4, 4);
			if (abs(cos(tmpalp)) < 0.0000001)
			{
				C(1, 1) = 0;
				C(1, 2) = -cos(tmpalp);
				C(2, 1) = cos(tmpalp);
				C(2, 2) = 0;
			}
			else if (abs(sin(tmpalp)) < 0.0000001)
			{
				C(1, 1) = -sin(tmpalp);
				C(1, 2) = 0;
				C(2, 1) = 0;
				C(2, 2) = -sin(tmpalp);
			}
			else
			{
				C(1, 1) = -sin(tmpalp);
				C(1, 2) = -cos(tmpalp);
				C(2, 1) = cos(tmpalp);
				C(2, 2) = -sin(tmpalp);
			}
			MatdAlpha.push_back(C);

			Matrix4d D = MatrixXd::Zero(4, 4);
			if (abs(cos(tmpthe)) < 0.0000001)
			{
				D(0, 0) = 0;
				D(0, 1) = -cos(tmpthe);
				D(1, 0) = cos(tmpthe);
				D(1, 1) = 0;
			}
			else if (abs(cos(tmpthe)) < 0.0000001)
			{
				D(0, 0) = -sin(tmpthe);
				D(0, 1) = 0;
				D(1, 0) = 0;
				D(1, 1) = -sin(tmpthe);
			}
			else
			{
				D(0, 0) = -sin(tmpthe);
				D(0, 1) = -cos(tmpthe);
				D(1, 0) = cos(tmpthe);
				D(1, 1) = -sin(tmpthe);
			}
			MatdTheta.push_back(D);

			Matrix4d A = MatrixXd::Zero(4, 4);
			A(0, 3) = 1;
			MatdA.push_back(A);

			Matrix4d B = MatrixXd::Zero(4, 4);
			B(2, 3) = 1;
			MatdD.push_back(B);
		}

	}
	
}

int DH::UpdateDiffJoint(vector<double>& theta, vector<vector<double>>& std)
{
	if (!ErrorModelFinishedBtn || theta.size() != 6 || std[0].size() != 6 || std[1].size() != 6 || std[2].size() != 6 || std[3].size() != 6)
		assert(0);

	for (int i = 0; i < std[0].size(); i++)
	{
		Matrix4d A = MatrixXd::Zero(4, 4);
		A(0, 3) = 1;
		MatdA[i] = A;

		Matrix4d B = MatrixXd::Zero(4, 4);
		B(2, 3) = 1;
		MatdD[i] = B;

		double tmpalp = twist_alpha[i] + std[2][i];
		Matrix4d C = MatrixXd::Zero(4, 4);
		if (abs(sin(tmpalp)) < 0.0000001)
		{
			C(1, 1) = 0;
			C(1, 2) = -cos(tmpalp);
			C(2, 1) = cos(tmpalp);
			C(2, 2) = 0;
		}
		else if (abs(cos(tmpalp)) < 0.0000001)
		{
			C(1, 1) = -sin(tmpalp);
			C(1, 2) = 0;
			C(2, 1) = 0;
			C(2, 2) = -sin(tmpalp);
		}
		else
		{
			C(1, 1) = -sin(tmpalp);
			C(1, 2) = -cos(tmpalp);
			C(2, 1) = cos(tmpalp);
			C(2, 2) = -sin(tmpalp);
		}
		MatdAlpha[i] = C;

		double tmpthe = angle_theta[i] + std[3][i] + theta[i] / 180 * pi;
		Matrix4d D = MatrixXd::Zero(4, 4);

		if (abs(sin(tmpthe)) < 0.0000001)
		{
			D(0, 0) = 0;
			D(0, 1) = -cos(tmpthe);
			D(1, 0) = cos(tmpthe);
			D(1, 1) = 0;
		}
		else if (abs(cos(tmpthe)) < 0.0000001)
		{
			D(0, 0) = -sin(tmpthe);
			D(0, 1) = 0;
			D(1, 0) = 0;
			D(1, 1) = -sin(tmpthe);
		}
		else
		{
			D(0, 0) = -sin(tmpthe);
			D(0, 1) = -cos(tmpthe);
			D(1, 0) = cos(tmpthe);
			D(1, 1) = -sin(tmpthe);
		}
		MatdTheta[i] = D;

	}
	return 0;
}



//////////////////////////////////////
// Calculation Function
/////////////////////////////////////

// Calculate End Pose by DH model
int DH::UpdateJoint(vector<double>& theta)
{
	if (theta.size() != 6)
		return -1;

	for (int i = 0; i < theta.size(); i++)
	{
		double angle = angle_theta[i] + theta[i] / 180.0 * pi;

		Matrix4d D = MatrixXd::Identity(4, 4);
		if (abs(cos(angle)) < 0.0000001)
		{
			D(0, 0) = 0;
			D(0, 1) = -sin(angle);
			D(1, 0) = sin(angle);
			D(1, 1) = 0;
		}
		else if (abs(sin(angle)) < 0.0000001)
		{
			D(0, 0) = cos(angle);
			D(0, 1) = 0;
			D(1, 0) = 0;
			D(1, 1) = cos(angle);
		}
		else
		{
			D(0, 0) = cos(angle);
			D(0, 1) = -sin(angle);
			D(1, 0) = sin(angle);
			D(1, 1) = cos(angle);
		}
		MatTheta[i] = D;

		Matrix4d E = MatAlpha[i] * MatA[i] * MatTheta[i] * MatD[i];
		MatJ[i] = E;
	}
}

Matrix4d DH::CalEndPose() {

	Eigen::Matrix4d P;

	P = MatJ[0] * MatJ[1] * MatJ[2] * MatJ[3] * MatJ[4] * MatJ[5] ;
	
	return P;
}

vector<Matrix4d> DH::CalEndPose(vector<vector<double>>& Theta)
{
	vector<Matrix4d> EndPose;

	for (int i = 0; i < Theta.size(); i++)
	{
		UpdateJoint(Theta[i]);
		Eigen::Matrix4d P1 = CalEndPose();
		EndPose.push_back(P1);
	}

	return EndPose;
}



// Calculate End Pose by MDH model
int DH::UpdateJointWithError(vector<double>& theta, vector<vector<double>>& std)
{
	if (!ErrorModelFinishedBtn || theta.size() != 6 || std[0].size() != 6 || std[1].size() != 6 || std[2].size() != 6 || std[3].size() != 6)
		assert(0);
		
	for (int i = 0; i < std[0].size(); i++)
	{
		Matrix4d A = MatrixXd::Identity(4, 4);
		A(0, 3) = length_a[i] + std[0][i];
		MateA[i] = A;

		Matrix4d B = MatrixXd::Identity(4, 4);
		B(2, 3) = offset_d[i] + std[1][i];
		MateD[i] = B;

		double tmpalp = twist_alpha[i] + std[2][i];
		Matrix4d C = MatrixXd::Identity(4, 4);
		if (abs(cos(tmpalp)) < 0.0000001)
		{
			C(1, 1) = 0;
			C(1, 2) = -sin(tmpalp);
			C(2, 1) = sin(tmpalp);
			C(2, 2) = 0;
		}
		else if (abs(sin(tmpalp)) < 0.0000001)
		{
			C(1, 1) = cos(tmpalp);
			C(1, 2) = 0;
			C(2, 1) = 0;
			C(2, 2) = cos(tmpalp);
		}
		else
		{
			C(1, 1) = cos(tmpalp);
			C(1, 2) = -sin(tmpalp);
			C(2, 1) = sin(tmpalp);
			C(2, 2) = cos(tmpalp);
		}
		MateAlpha[i] = C;

		double tmpthe = angle_theta[i] + std[3][i] + theta[i] / 180 * pi;
		Matrix4d D = MatrixXd::Identity(4, 4);

		if (abs(cos(tmpthe)) < 0.0000001)
		{
			D(0, 0) = 0;
			D(0, 1) = -sin(tmpthe);
			D(1, 0) = sin(tmpthe);
			D(1, 1) = 0;
		}
		else if (abs(sin(tmpthe)) < 0.0000001)
		{
			D(0, 0) = cos(tmpthe);
			D(0, 1) = 0;
			D(1, 0) = 0;
			D(1, 1) = cos(tmpthe);
		}
		else
		{
			D(0, 0) = cos(tmpthe);
			D(0, 1) = -sin(tmpthe);
			D(1, 0) = sin(tmpthe);
			D(1, 1) = cos(tmpthe);
		}
		MateTheta[i] = D;

		Matrix4d E = MateAlpha[i] * MateA[i] * MateTheta[i] * MateD[i];
		MateJ[i] = E;
	}
	return 0;
}

Matrix4d DH::CalEndPoseWithError()
{
	Matrix4d eP;

	eP = MateJ[0] * MateJ[1] * MateJ[2] * MateJ[3] * MateJ[4] * MateJ[5];

	return eP;
}

vector<Matrix4d> DH::CalEndPoseWithError(vector<vector<double>>& Theta, vector<vector<double>>& std) {
	
	vector<Matrix4d> EndPose;
	if (!DH::ErrorModelFinishedBtn)
		CreateJointWithError(std);

	for (int i = 0; i < Theta.size(); i++)
	{
		UpdateJointWithError(Theta[i], std);
		Eigen::Matrix4d P1 = CalEndPoseWithError();
		EndPose.push_back(P1);
	}

	return EndPose;
}



vector<Vector4d> DH::CalTargetPose(vector<Matrix4d>& EndPose, Vector4d b) {
	vector<Vector4d> TargetPose;

	for (int i = 0; i < EndPose.size(); i++)
	{
		Vector4d tmp = EndPose[i] * b;
		TargetPose.push_back(tmp);
	}

	return TargetPose;
}

vector<Vector3d> DH::CalDiff_Tar2Laser(vector<Vector4d> TargetPose,vector<Vector4d> Laser, Matrix4d Trans_Tar2Laser) {
	
	
	if (Laser.size() != TargetPose.size())
		assert(1);
	diff_all_T2L.resize(Laser.size());
	
	for (int i = 0; i < Laser.size(); i++)
	{
		TargetPose[i] = Trans_Tar2Laser * TargetPose[i];
		
		Vector3d TargetPose3;
		TargetPose3[0] = TargetPose[i][0] / TargetPose[i][3];
		TargetPose3[1] = TargetPose[i][1] / TargetPose[i][3];
		TargetPose3[2] = TargetPose[i][2] / TargetPose[i][3];

		Vector3d Laser3;
		Laser3[0] = Laser[i][0] / Laser[i][3];
		Laser3[1] = Laser[i][1] / Laser[i][3];
		Laser3[2] = Laser[i][2] / Laser[i][3];

		Vector3d tmp = Laser3 - TargetPose3;
		diff_all_T2L[i] = tmp;
	}

	return diff_all_T2L;
}

vector<Vector3d> DH::CalDiff_Laser2Tar(vector<Vector4d> Laser, vector<Vector4d> TargetPose, Matrix4d Trans_Tar2Laser) {


	if (Laser.size() != TargetPose.size())
		assert(1);
	diff_all_L2T.resize(Laser.size());

	for (int i = 0; i < Laser.size(); i++)
	{
		TargetPose[i] = Trans_Tar2Laser * TargetPose[i];

		Vector3d TargetPose3;
		TargetPose3[0] = TargetPose[i][0] / TargetPose[i][3];
		TargetPose3[1] = TargetPose[i][1] / TargetPose[i][3];
		TargetPose3[2] = TargetPose[i][2] / TargetPose[i][3];

		Vector3d Laser3;
		Laser3[0] = Laser[i][0] / Laser[i][3];
		Laser3[1] = Laser[i][1] / Laser[i][3];
		Laser3[2] = Laser[i][2] / Laser[i][3];

		Vector3d tmp = TargetPose3 - Laser3;
		diff_all_L2T[i] = tmp;
	}

	return diff_all_L2T;
}




void DH::SumDiff(vector<Vector3d> Diff) {
	DH::Sum_Diff = 0;

	for (int i = 0; i < Diff.size(); i++)
	{
		Sum_Diff += Diff[i].transpose() * Diff[i];
	}
}

MatrixXd DH::TransDiff2Fp(vector<Vector3d>& TmpDiff) {

	int num = 3 * TmpDiff.size();
	MatrixXd Fp;
	Fp.resize(num,1);

	for (int i = 0; i < TmpDiff.size(); i++)
	{
		int j = i * 3;
		Fp.block<3, 1>(j, 0) = TmpDiff[i];
	}

	return Fp;

}

Matrix4d DH::dAlp(int i) {
	return MatdAlpha[i] * MateA[i] * MateTheta[i] * MateD[i];
}
Matrix4d DH::dA(int i) {
	return MateAlpha[i] * MatdA[i] * MateTheta[i] * MateD[i];
}
Matrix4d DH::dThe(int i) {
	return MateAlpha[i] * MateA[i] * MatdTheta[i] * MateD[i];
}
Matrix4d DH::dD(int i) {
	return MateAlpha[i] * MateA[i] * MateTheta[i] * MatdD[i];
}


int DH::MDH_Solver(double Sum_Delta, int MaxIter, double MinError, double MinDelta,
	vector<vector<double>> Theta, vector<Vector4d> Laser, Vector4d b1, Matrix4d Trans_Tar2Laser) 
{
	int iter = 0;
	vector<vector<double>> std = { { 0,0,0,0,0,0 },{ 0,0,0,0,0,0 },{ 0,0,0,0,0,0 },{ 0,0,0,0,0,0 } };

	CreateDiffJoint(std);
	vector<Matrix4d> EndPose = CalEndPoseWithError(Theta, std);

	vector<Vector4d> TargetPose = CalTargetPose(EndPose, b1);
	vector<Vector3d> TmpDiff = CalDiff_Laser2Tar( TargetPose, Laser, Trans_Tar2Laser);

	MatrixXd Fp = TransDiff2Fp(TmpDiff);
	
	SumDiff(TmpDiff);
	int num = Theta.size() * 3;
	while (Sum_Diff > MinError && Sum_Delta > MinDelta && iter <= MaxIter) {
		
		MatrixXd Jacobi;
		Jacobi.resize(num, 23);
		for (int i = 0; i < Theta.size(); i++)
		{
			UpdateDiffJoint(Theta[i], std);
			UpdateJointWithError(Theta[i], std);
			Matrix<double, 4, 6> dTdA, dTdd, dTdAlp, dTdThe;
			
			dTdA.col(0) = Trans_Tar2Laser * dA(0) * MateJ[1] * MateJ[2] * MateJ[3] * MateJ[4] * MateJ[5] * b1;
			dTdd.col(0) = Trans_Tar2Laser * dD(0) * MateJ[1] * MateJ[2] * MateJ[3] * MateJ[4] * MateJ[5] * b1;
			dTdAlp.col(0) = Trans_Tar2Laser * dAlp(0) * MateJ[1] * MateJ[2] * MateJ[3] * MateJ[4] * MateJ[5] * b1;
			dTdThe.col(0) = Trans_Tar2Laser * dThe(0) * MateJ[1] * MateJ[2] * MateJ[3] * MateJ[4] * MateJ[5] * b1;

			dTdA.col(1) = Trans_Tar2Laser * MateJ[0] * dA(1) * MateJ[2] * MateJ[3] * MateJ[4] * MateJ[5] * b1;
			dTdd.col(1) = Trans_Tar2Laser *  MateJ[0] * dD(1) * MateJ[2] * MateJ[3] * MateJ[4] * MateJ[5] * b1;
			dTdAlp.col(1) = Trans_Tar2Laser * MateJ[0] * dAlp(1) *  MateJ[2] * MateJ[3] * MateJ[4] * MateJ[5] * b1;
			dTdThe.col(1) = Trans_Tar2Laser * MateJ[0] * dThe(1) *  MateJ[2] * MateJ[3] * MateJ[4] * MateJ[5] * b1;

			dTdA.col(2) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * dA(2) *  MateJ[3] * MateJ[4] * MateJ[5] * b1;
			//dTdd.col(2) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * dD(2) *  MateJ[3] * MateJ[4] * MateJ[5] * b1;
			dTdAlp.col(2) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * dAlp(2) *  MateJ[3] * MateJ[4] * MateJ[5] * b1;
			dTdThe.col(2) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * dThe(2) *  MateJ[3] * MateJ[4] * MateJ[5] * b1;

			dTdA.col(3) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * MateJ[2] * dA(3) * MateJ[4] * MateJ[5] * b1;
			dTdd.col(2) = Trans_Tar2Laser * MateJ[0] * MateJ[1] *   MateJ[2] * dD(3) * MateJ[4] * MateJ[5] * b1;
			dTdAlp.col(3) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * MateJ[2] * dAlp(3) * MateJ[4] * MateJ[5] * b1;
			dTdThe.col(3) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * MateJ[2] * dThe(3) * MateJ[4] * MateJ[5] * b1;

			dTdA.col(4) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * MateJ[2] * MateJ[3] * dA(4) * MateJ[5] * b1;
			dTdd.col(3) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * MateJ[2] * MateJ[3] * dD(4) * MateJ[5] * b1;
			dTdAlp.col(4) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * MateJ[2] * MateJ[3] * dAlp(4) * MateJ[5] * b1;
			dTdThe.col(4) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * MateJ[2] * MateJ[3] * dThe(4) * MateJ[5] * b1;

			dTdA.col(5) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * MateJ[2] * MateJ[3] *  MateJ[4] * dA(5) *b1;
			dTdd.col(4) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * MateJ[2] * MateJ[3] *  MateJ[4] * dD(5) *b1;
			dTdAlp.col(5) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * MateJ[2] * MateJ[3] * MateJ[4] * dAlp(5) * b1;
			dTdThe.col(5) = Trans_Tar2Laser * MateJ[0] * MateJ[1] * MateJ[2] * MateJ[3] * MateJ[4] * dThe(5) * b1;

			Matrix<double, 3, 23> JcbTmp;
			JcbTmp.block<3, 6>(0, 0) = dTdA.block<3, 6>(0, 0);
			JcbTmp.block<3, 5>(0, 6) = dTdd.block<3, 5>(0, 0);
			JcbTmp.block<3, 6>(0, 11) = dTdAlp.block<3, 6>(0, 0);
			JcbTmp.block<3, 6>(0, 17) = dTdThe.block<3, 6>(0, 0);

					
			

			int j = i * 3;
			Jacobi.block<3, 23>(j, 0) = JcbTmp;

		}

		MatrixXd PLInfo = Jacobi.transpose() * Jacobi;
		MatrixXd JcbB = -Jacobi.transpose() * Fp;
		VectorXd Delta = PLInfo.inverse() * JcbB;
		Sum_Delta = Delta.transpose() * Delta;

		// a
		std[0][0] += Delta[0];
		std[0][1] += Delta[1];
		std[0][2] += Delta[2];
		std[0][3] += Delta[3];
		std[0][4] += Delta[4];
		std[0][5] += Delta[5];
		// d
		std[1][0] += Delta[6];
		std[1][1] += Delta[7];
		std[1][2] += 0;
		std[1][3] += Delta[8];
		std[1][4] += Delta[9];
		std[1][5] += Delta[10];
		// alp
		std[2][0] += Delta[11];
		std[2][1] += Delta[12];
		std[2][2] += Delta[13];
		std[2][3] += Delta[14];
		std[2][4] += Delta[15];
		std[2][5] += Delta[16];
		// the
		std[3][0] += Delta[17];
		std[3][1] += Delta[18];
		std[3][2] += Delta[19];
		std[3][3] += Delta[20];
		std[3][4] += Delta[21];
		std[3][5] += Delta[22];

		EndPose = CalEndPoseWithError(Theta, std);

		TargetPose = CalTargetPose(EndPose, b1);
		TmpDiff = CalDiff_Laser2Tar(Laser, TargetPose, Trans_Tar2Laser);

		Fp = TransDiff2Fp(TmpDiff);

		SumDiff(TmpDiff);
		
		++iter;
		cout<<"Iterations "<< iter<<" Error "<< Sum_Diff << endl;
	}

	return 0;
}


///////////////////////////////////////
// Private Member Read from TXT 
int DH::ReadDHpara(FILE* fp)
{
	int DHNo = 0, iter = 0, ch;
	

	while (!feof(fp)) {
		if ((ch = fgetc(fp)) == '#') { /* skip comments */
			SKIP_LINE(fp);
			continue;
		}

		if (feof(fp))
		{
			break;
		}
		double tmp[6] = {0,0,0,0,0,0};

		int nret = fscanf(fp, "%lf %lf %lf %lf %lf %lf", tmp, tmp + 1, tmp + 2, tmp + 3, tmp + 4, tmp + 5);
		if (nret == -1)
			break;
		if (nret == 6)
		{
			++iter;

			if (iter == 1) 	for (int i = 0; i < 6; i++)
				length_a[i] = tmp[i];
			if (iter == 2) 	for (int i = 0; i < 6; i++)
				twist_alpha[i] = tmp[i] / 180.0 * pi;
			if (iter == 3) 	for (int i = 0; i < 6; i++)
				offset_d[i] = tmp[i];
			if (iter == 4) 	for (int i = 0; i < 6; i++)
				angle_theta[i] = tmp[i] / 180.0 * pi;

			SKIP_LINE(fp);
		}
		

		DHNo++;
		}

		
		return DHNo;
	}

void DH::CreateInitJoint() {

	for (int i = 0; i < 6; i++)
	{
		Matrix4d A = MatrixXd::Identity(4, 4);
		A(0, 3) = length_a[i];
		MatA.push_back(A);

		Matrix4d B = MatrixXd::Identity(4, 4);
		B(2, 3) = offset_d[i];
		MatD.push_back(B);

		Matrix4d C = MatrixXd::Identity(4, 4);
		if (abs(cos(twist_alpha[i])) < 0.0000001)
		{
			C(1, 1) = 0;
			C(1, 2) = -sin(twist_alpha[i]);
			C(2, 1) = sin(twist_alpha[i]);
			C(2, 2) = 0;
		}
		else if (abs(sin(twist_alpha[i])) < 0.0000001)
		{
			C(1, 1) = cos(twist_alpha[i]);
			C(1, 2) = 0;
			C(2, 1) = 0;
			C(2, 2) = cos(twist_alpha[i]);
		}
		else
		{
			C(1, 1) = cos(twist_alpha[i]);
			C(1, 2) = -sin(twist_alpha[i]);
			C(2, 1) = sin(twist_alpha[i]);
			C(2, 2) = cos(twist_alpha[i]);
		}
		MatAlpha.push_back(C);

		Matrix4d D = MatrixXd::Identity(4, 4);
		if (abs(cos(angle_theta[i])) < 0.0000001)
		{
			D(0, 0) = 0;
			D(0, 1) = -sin(angle_theta[i]);
			D(1, 0) = sin(angle_theta[i]);
			D(1, 1) = 0;
		}
		else if (abs(sin(angle_theta[i])) < 0.0000001)
		{
			D(0, 0) = cos(angle_theta[i]);
			D(0, 1) = 0;
			D(1, 0) = 0;
			D(1, 1) = cos(angle_theta[i]);
		}
		else
		{
			D(0, 0) = cos(angle_theta[i]);
			D(0, 1) = -sin(angle_theta[i]);
			D(1, 0) = sin(angle_theta[i]);
			D(1, 1) = cos(angle_theta[i]);
		}
		MatTheta.push_back(D);

		Matrix4d E = MatAlpha[i] * MatA[i] * MatTheta[i] * MatD[i];
		MatJ.push_back(E);
	}
	
}


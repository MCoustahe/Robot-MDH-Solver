// Robot.cpp : Defines the entry point for the console application.
//
//

#include "stdafx.h"
#include "string"
#include "math.h"
#include <cmath>

#include "DHClass.h"


using namespace std;
using namespace Eigen;


void PrintMatrix4d(Eigen::Matrix4d P)
{
	
	for (int j = 0; j < 4; j++) {
		for (int k = 0; k < 4; k++) cout << P(j, k) << " ";
		cout << endl;
	}
	
}

vector<vector<double>> ReadTheta(FILE* fp)
{
	vector<vector<double>> Theta ;
	
	while (!feof(fp))
	{
		
		double theta1[6] = {0,0,0,0,0,0};
		int nret = fscanf(fp, "%lf %lf %lf %lf %lf %lf", theta1, theta1 + 1, theta1 + 2, theta1 + 3, theta1 + 4, theta1 + 5);
		if (nret == -1)
			break;

		vector<double> tmp;
		if (nret == 6)
		{
			for (int j = 0; j < 6; j++)
				tmp.push_back(theta1[j]);
			
		}
		Theta.push_back(tmp);

	}

	return Theta;
}

vector<Vector4d> ReadLaser(FILE* fp){

	vector<Vector4d> Laser;

	while (!feof(fp))
	{
		int num = 0;
		long double L[3] = { 0,0,0 };
		int nret = fscanf(fp, "%d %lf %lf %lf", &num, L, L + 1, L + 2);
		if (nret == -1)
			break;

		Vector4d tmp;
		if (nret == 4)
		{
			tmp << L[0], L[1], L[2], 1;
		}
		Laser.push_back(tmp);

	}

	return Laser;
}


int _tmain(int argc, _TCHAR* argv[])
{
	DH robot(true);

	// Read Laser data and joint angles 
	char* szTheta = "../data/Theta.txt";
	char* szLaser = "../data/laser.txt";
	FILE* fpTheta = fopen(szTheta, "r");
	FILE* fpLaser = fopen(szLaser, "r");
	vector<vector<double>> Theta = ReadTheta(fpTheta);
	vector<Vector4d> Laser = ReadLaser(fpLaser);

	// Transformation between the end and target
	Vector4d b1 = {25.7385321859648,83.0153073423998,30.4584412221448,1};
	
	// Calculate by DH model
	vector<Matrix4d> EndPose = robot.CalEndPose(Theta);

	vector<Vector4d> TargetPose = robot.CalTargetPose(EndPose, b1);

	// Transformation between Laser coordinate system and robot coordinate system
	Matrix4d Trans_Tar2Laser;
	Trans_Tar2Laser << -0.987085968911664, - 0.160069976387521, 0.00623639615629754,	4331.12233008149,
						0.160090725398695, - 0.987097805383778, 0.00298031002551424, - 4832.1141187539,
						0.0056788748039735, 0.00394021139372697,0.999976112272255, - 1205.20280275572,
						0,                      0,               0,                  1 ;

	robot.MDH_Solver(22, 50, 1e-8, 1e-10, Theta, Laser, b1, Trans_Tar2Laser);

	while (1);

	return 0;
}


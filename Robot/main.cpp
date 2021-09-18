#include "matrix.h"
#include "DHClass.h"
#include "fstream"
#include <string>
using namespace std;


double theta0[6] = { 0,0,0,0,0,0 };
double theta1[6] = { 29.34,
					-6.43,
					28.35,
					32.69,
					-16.31,
					-22.74 };
double theta2[6] = { 1.38,
					18.87,
					-3.74,
					5.39,
					57.95,
					0.41 };
double theta3[6] = { 0.53,
					38.28,
					-9.06,
					7.45,
					13.61,
					-2.48 }; 
double theta4[6] = { -2.07,
					1.4,
					43.56,
					-24.48,
					10.31,
					-67.22 };
double theta5[6] = { -20.69,
					-15.74,
					44.15,
					-25.01,
					32.23,
					-67.22 };
double theta6[6] = { -12.94,
					-35.93,
					39.25,
					-56.62,
					-14.30,
					-29.84 };
double theta7[6] = { -10.1,
					10.25,
					63.24,
					-9.27,
					-78.26,
					-149.45 };
double theta8[6] = { 19.33,	
					31.48,
					28.64,
					22.17,
					-40.19,
					-43.30 };


double loc1[3] = { 1343.69,720.35,1629.08 };
double loc2[3] = { 1896.38,61.6,1548.81 };
double loc3[3] = { 2258.43,27.02,1144.15 };
double loc4[3] = { 1411.94,-65.75,1076.82 };
double loc5[3] = { 1123.25,-472.36,1329.79 };
double loc6[3] = { 990.82,-185.29,1840.30 };
double loc7[3] = { 1219.80,-185.32,864.29 };
double loc8[3] = { 1742.85,559.75,775.48 };

double q1[4] = { 0.636347,-0.0846056,0.749928,0.159726 };
double q2[4] = { 0.146292,0.0104785,0.988146,0.0453693 };
double q3[4] = { 0.399861,0.0277232,0.91558,0.0324778 };
double q4[4] = { 0.184033,-0.660405,0.684891,-0.246823 };
double q5[4] = { 0.0683984,-0.508897,0.809192,-0.285577 };
double q6[4] = { 0.530294,-0.460926,0.497822,-0.508438 };
double q7[4] = { 0.170984,-0.656947,0.169389,-0.714487 };
double q8[4] = { 0.547518,-0.254032,0.785983,-0.133872 };

double MJ1[16] = Z, MJ2[16] = Z, MJ3[16] = Z, MJ4[16] = Z, MJ5[16] = Z, MJ6[16] = Z, MJ7[16] = Z, MJ8[16] = Z;

void CreateJMatrix_test()
{
	CreateJointMatrix(q1, loc1, MJ1);
	CreateJointMatrix(q2, loc2, MJ2);
	CreateJointMatrix(q3, loc3, MJ3);
	CreateJointMatrix(q4, loc4, MJ4);
	CreateJointMatrix(q5, loc5, MJ5);
	CreateJointMatrix(q6, loc6, MJ6);
	CreateJointMatrix(q7, loc7, MJ7);
	CreateJointMatrix(q8, loc8, MJ8);

	if (PRINT_ALL||PRINT_MJ_ALL)
	{
		PrintMatrix_OneDim_double(MJ1);
		PrintMatrix_OneDim_double(MJ2);
		PrintMatrix_OneDim_double(MJ3);
		PrintMatrix_OneDim_double(MJ4);
		PrintMatrix_OneDim_double(MJ5);
		PrintMatrix_OneDim_double(MJ6);
		PrintMatrix_OneDim_double(MJ7);
		PrintMatrix_OneDim_double(MJ8);
	}

}

void FK_test()
{
	DH ABB6700(length_a_6700_200_260_4, twist_alpha_6700_200_260_4,
			   offset_d_6700_200_260_4, angle_theta_6700_200_260_4);
	double m1[16] = Z, m2[16] = Z, m3[16] = Z, m4[16] = Z, m5[16] = Z, m6[16] = Z, m7[16] = Z, m8[16] = Z;

	//ABB6700.Forward_Kinematics(theta0, m1);
	ABB6700.Forward_Kinematics(theta1, m1);
	ABB6700.Forward_Kinematics(theta2, m2);
	ABB6700.Forward_Kinematics(theta3, m3);
	ABB6700.Forward_Kinematics(theta4, m4);
	ABB6700.Forward_Kinematics(theta5, m5);
	ABB6700.Forward_Kinematics(theta6, m6);
	ABB6700.Forward_Kinematics(theta7, m7);
	ABB6700.Forward_Kinematics(theta8, m8);

	if (PRINT_ALL||PRINT_FK_ALL)
	{
		PrintMatrix_OneDim_double(m1);
		PrintMatrix_OneDim_double(m2);
		PrintMatrix_OneDim_double(m3);
		PrintMatrix_OneDim_double(m4);
		PrintMatrix_OneDim_double(m5);
		PrintMatrix_OneDim_double(m6);
		PrintMatrix_OneDim_double(m7);
		PrintMatrix_OneDim_double(m8);
	}

}

void IK_test(FILE* fpTheta)
{
	DH ABB6700(length_a_6700_200_260_4, twist_alpha_6700_200_260_4,
			   offset_d_6700_200_260_4, angle_theta_6700_200_260_4);


	double  thetainv1[6] = { 0,0,0,0,0,0 },
			thetainv2[6] = { 0,0,0,0,0,0 }, 
			thetainv3[6] = { 0,0,0,0,0,0 }, 
			thetainv4[6] = { 0,0,0,0,0,0 }, 
			thetainv5[6] = { 0,0,0,0,0,0 }, 
			thetainv6[6] = { 0,0,0,0,0,0 }, 
			thetainv7[6] = { 0,0,0,0,0,0 }, 
			thetainv8[6] = { 0,0,0,0,0,0 };

	ABB6700.Inverse_Kinematics(MJ1, thetainv1);
	ABB6700.Inverse_Kinematics(MJ2, thetainv2);
	ABB6700.Inverse_Kinematics(MJ3, thetainv3);
	ABB6700.Inverse_Kinematics(MJ4, thetainv4);
	ABB6700.Inverse_Kinematics(MJ5, thetainv5);
	ABB6700.Inverse_Kinematics(MJ6, thetainv6);
	ABB6700.Inverse_Kinematics(MJ7, thetainv7);
	ABB6700.Inverse_Kinematics(MJ8, thetainv8);

	fprintf(fpTheta, "%lf %lf %lf %lf %lf %lf\n", thetainv1[0], thetainv1[1], thetainv1[2], thetainv1[3], thetainv1[4], thetainv1[5]);
	fprintf(fpTheta, "%lf %lf %lf %lf %lf %lf\n", thetainv2[0], thetainv2[1], thetainv2[2], thetainv2[3], thetainv2[4], thetainv2[5]);
	fprintf(fpTheta, "%lf %lf %lf %lf %lf %lf\n", thetainv3[0], thetainv3[1], thetainv3[2], thetainv3[3], thetainv3[4], thetainv3[5]);
	fprintf(fpTheta, "%lf %lf %lf %lf %lf %lf\n", thetainv4[0], thetainv4[1], thetainv4[2], thetainv4[3], thetainv4[4], thetainv4[5]);
	fprintf(fpTheta, "%lf %lf %lf %lf %lf %lf\n", thetainv5[0], thetainv5[1], thetainv5[2], thetainv5[3], thetainv5[4], thetainv5[5]);
	fprintf(fpTheta, "%lf %lf %lf %lf %lf %lf\n", thetainv6[0], thetainv6[1], thetainv6[2], thetainv6[3], thetainv6[4], thetainv6[5]);
	fprintf(fpTheta, "%lf %lf %lf %lf %lf %lf\n", thetainv7[0], thetainv7[1], thetainv7[2], thetainv7[3], thetainv7[4], thetainv7[5]);
	fprintf(fpTheta, "%lf %lf %lf %lf %lf %lf\n", thetainv8[0], thetainv8[1], thetainv8[2], thetainv8[3], thetainv8[4], thetainv8[5]);

	
	if (PRINT_ALL || PRINT_IK_ALL)
	{
		PrintMatrix_OneDim_double6(thetainv1);
		PrintMatrix_OneDim_double6(thetainv2);
		PrintMatrix_OneDim_double6(thetainv3);
		PrintMatrix_OneDim_double6(thetainv4);
		PrintMatrix_OneDim_double6(thetainv5);
		PrintMatrix_OneDim_double6(thetainv6);
		PrintMatrix_OneDim_double6(thetainv7);
		PrintMatrix_OneDim_double6(thetainv8);
	}
	
	fclose(fpTheta);
}

int main()
{
	const char* szTheta = "../data/theta.txt";
	FILE* fpTheta = fopen(szTheta, "w");

	//ABB6700.getlength_a();
	//ABB6700.gettwist_alpha();
	//ABB6700.getoffset_d();
	//ABB6700.getangle_theta();
	//vector<double> ve(16, 0);
	//vector<double> ve1(3, 0);
	//ve1.push_back(1);
	//ve1.push_back(0);

	//PrintMatrix_OneDim_double1(ve);
	//PrintMatrix_OneDim_double1(ve1);

	cout << "DH test : " << endl;
	cout << endl;


	cout << " CreateEnd-EFFECTOR : " << endl;
	cout << endl;
	CreateJMatrix_test();
	cout << " Forward Kinematic : " << endl;
	cout << endl;
	FK_test();
	cout << " Inverse Kinematic : " << endl;
	cout << endl;
	//IK_test(fpTheta);


	//ABB6700.Inverse_Kinematics(m, theta);
	//cout << " theta[0] : " << theta[0] << endl;
	//cout << " theta[1] : " << theta[1] << endl;
	//cout << " theta[2] : " << theta[2] << endl;
	//cout << " theta[3] : " << theta[3] << endl;
	//cout << " theta[4] : " << theta[4] << endl;
	//cout << " theta[5] : " << theta[5] << endl;

	//Forward_Kinematics(theta1);
	//Forward_Kinematics(theta2);
	//Forward_Kinematics(theta3);
	//Forward_Kinematics(theta4);
	//Forward_Kinematics(theta5);
	//Forward_Kinematics(theta6);
	//Forward_Kinematics(theta7);
	//Forward_Kinematics(theta8);

	
	//CreateQuaternionMatrix(0, 0.5, 0.86603, 0, q);

	//position2
	//CreateQuaternionMatrix(0.146292, 0.0104785, 0.988145, 0.0453693, q);
	////position3
	//CreateQuaternionMatrix(0.399861, 0.0277232, 0.91558, 0.03324778, q);
	//PrintMatrix_OneDim_double(q);
	
	/*
	double* A = new double[N * N]();
	srand((unsigned)time(0));
	for (int i = 0; i < N; i++)
	{
	    for (int j = 0; j < N; j++)
	    {
	        A[i * N + j] = rand() % 100 * 0.01;
	    }
	}
	
	
	double* E_test = new double[N * N]();
	double* invOfA = new double[N * N]();
	invOfA = LUP_solve_inverse(A);
	
	E_test = mul(A, invOfA);    //验证精确度
	
	cout << "矩阵A:" << endl;
	for (int i = 0; i < N; i++)
	{
	    for (int j = 0; j < N; j++)
	    {
	        cout << A[i * N + j] << " ";
	    }
	    cout << endl;
	}
	
	cout << "inv_A:" << endl;
	for (int i = 0; i < N; i++)
	{
	    for (int j = 0; j < N; j++)
	    {
	        cout << invOfA[i * N + j] << " ";
	    }
	    cout << endl;
	}
	
	cout << "E_test:" << endl;
	for (int i = 0; i < N; i++)
	{
	    for (int j = 0; j < N; j++)
	    {
	        cout << E_test[i * N + j] << " ";
	    }
	    cout << endl;
	}

	delete A;
	delete E_test;
	delete invOfA;
	*/
	return 0;
}
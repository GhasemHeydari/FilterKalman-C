#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// This project has developed for UKF algorithm for estimation of a orbit.
// In this program I had tried to code many function that my supplier wanted.

double InputRealPositionInECI[3];
double InputQuatOfSatECIToBody[4];
double mu;
double EstimateOfR[3], EstimateOfV[3];
double lambda, c, gamma, alpha, beta, ki;
double omega, ZHat, Pzz, MeasurmentNoise;
double A[6][6];
double Y1[6][13], Chi[6][13], ChiMinus[6][13], ChiPlus[6][13], YY[6][13];
double Y3[6][13];
double PMinus[6][6], PMinusP[6][6];
double DiagonalOmegaC[13][13];
double ProccessNoise[6][6], xHatMinus[6][13];
double OmegaS[13], OmegaC[13];
double xHat[6], xHatM[6], xHatP[6], GammaK[13], Pxz[6], k[6];
double StateEstimate[6], Y2[13];
double state[6];
double InitialError[6];
double InputCurrentMeasurment;
double MultiNormMagnetInNED[13];
int n;
int m;
int TimeMission;
int firstTime;
///
void Initialization();
void RungeKutta4(double, int, double[], double, double, double[]);
void Predict();
void Update(double[], double);
void CalcSigma(double state[], double covariance[6][6], double gamma, int n, double Chi_Output[6][13]);

int main()
{

	// Time of Run
	TimeMission = 2;

	if (firstTime == 0)
	{
		/// Initialization function for set requirement settings
		Initialization();
		firstTime = 1;
	}
	

	for (int i = 0; i < TimeMission; i++)
	{
		if (i == 0)
		{
			// Inputs of every time
			InputQuatOfSatECIToBody[0] = 0;
			InputQuatOfSatECIToBody[1] = 0;
			InputQuatOfSatECIToBody[2] = 0;
			InputQuatOfSatECIToBody[3] = 1;

			InputCurrentMeasurment = 0;
			
			// Data of IGRF(That generate in Update method)
			MultiNormMagnetInNED[0] = 3.2559253926989378E-05;
			MultiNormMagnetInNED[1] = 3.0875829361735096E-05;
			MultiNormMagnetInNED[2] = 3.2654778982637390E-05;
			MultiNormMagnetInNED[3] = 3.2424533842081306E-05;
			MultiNormMagnetInNED[4] = 3.2559253926989378E-05;
			MultiNormMagnetInNED[5] = 3.2559253926989378E-05;
			MultiNormMagnetInNED[6] = 3.2559253926989378E-05;
			MultiNormMagnetInNED[7] = 3.4370504410974555E-05;
			MultiNormMagnetInNED[8] = 3.2433278428463369E-05;
			MultiNormMagnetInNED[9] = 3.2686083868830905E-05;
			MultiNormMagnetInNED[10] = 3.2559253926989378E-05;
			MultiNormMagnetInNED[11] = 3.2559253926989378E-05;
			MultiNormMagnetInNED[12] = 3.2559253926989378E-05; 

		}
		else
		{
			// Inputs of every time
			InputQuatOfSatECIToBody[0] = 0;
			InputQuatOfSatECIToBody[1] = 0;
			InputQuatOfSatECIToBody[2] = 0;
			InputQuatOfSatECIToBody[3] = 1;

			InputCurrentMeasurment = 3.299266053633676E-05;

			// Data of IGRF(That generate in Update method)
			MultiNormMagnetInNED[0]  = 1.4110772735650189E-05;
			MultiNormMagnetInNED[1]  = 1.4040929425511113E-05;
			MultiNormMagnetInNED[2]  = 1.4137474082594834E-05;
			MultiNormMagnetInNED[3]  = 1.4067976414198402E-05;
			MultiNormMagnetInNED[4]  = 1.4110772735650189E-05;
			MultiNormMagnetInNED[5]  = 1.4110772735650189E-05;
			MultiNormMagnetInNED[6]  = 1.4110772735650189E-05;
			MultiNormMagnetInNED[7]  = 1.4180212459972238E-05;
			MultiNormMagnetInNED[8]  = 1.4077378179769822E-05;
			MultiNormMagnetInNED[9]  = 1.4152617613816043E-05;
			MultiNormMagnetInNED[10] = 1.4110772735650189E-05;
			MultiNormMagnetInNED[11] = 1.4110772735650189E-05;
			MultiNormMagnetInNED[12] = 1.4110772735650189E-05;
		}

		// InputQuatOfSatECIToBody data used in dynamic(RK4) in drag calculation.
		Predict(InputQuatOfSatECIToBody);
		// InputCurrentMeasurment data with noise and NormMagnetInNED Data without noise that calculated in Filter.
		Update(MultiNormMagnetInNED, InputCurrentMeasurment);

		printf("xhat: %f,%f,%f,%f,%f,%f \n", xHat[0], xHat[1], xHat[2], xHat[3], xHat[4], xHat[5]);
	}

	return 0;
}
void Initialization()
{
	state[0] = 6788137.0199246109;
	state[1] = 0;
	state[2] = 0;
	state[3] = 0;
	state[4] = -7432.0239552857747;
	state[5] = -1901.2041812034477;
	InitialError[0] = 20000;
	InitialError[1] = 20000;
	InitialError[2] = 20000;
	InitialError[3] = 50;
	InitialError[4] = 50;
	InitialError[5] = 50;
	alpha 	= 1;
	beta 	= 2;
	ki 		= 1;
	double InitialState[6];
	InitialState[0] = state[0];
	InitialState[1] = state[1];
	InitialState[2] = state[2];
	InitialState[3] = state[3];
	InitialState[4] = state[4];
	InitialState[5] = state[5];

	n = 6;
	double StateEstimate[6];
	VectorZeros(n, StateEstimate);
	Vector_Sum(n, InitialState, InitialError, StateEstimate);

	EqualVector(n, StateEstimate, xHat);

	double InitialErrorpow2[6];
	VectorZeros(n, InitialErrorpow2);
	InitialErrorpow2[0] = pow(InitialError[0], 2);
	InitialErrorpow2[1] = pow(InitialError[1], 2);
	InitialErrorpow2[2] = pow(InitialError[2], 2);
	InitialErrorpow2[3] = pow(InitialError[3], 2);
	InitialErrorpow2[4] = pow(InitialError[4], 2);
	InitialErrorpow2[5] = pow(InitialError[5], 2);
	VectorZeros(n, Pxz);
	VectorZeros(n, k);
	Matrix_Zero(n, n, PMinus);
	for (int i = 0; i < n; i++) {
		PMinus[i][i] = 4 * InitialErrorpow2[i];
	}
	//----------------------------------------------------------------
	m = 1;
	lambda = (pow(alpha, 2) * (n + ki)) - n;
	c = lambda + n;
	gamma = sqrt(c);

	VectorZeros(1 + (2 * n), OmegaS);
	OmegaS[0] = lambda / c;
	for (int j = 1; j < (2 * n) + 1; j++)
	{
		OmegaS[j] = 1 / (2 * c);
	}
	OmegaC[0] = (lambda / c) + (1 - pow(alpha, 2) + beta);
	for (int j = 1; j < (2 * n) + 1; j++)
	{
		OmegaC[j] = 1 / (2 * c);
	}
	Matrix_Zero(1 + (2 * n), 1 + (2 * n), DiagonalOmegaC);
	for (int i = 0; i < 1 + (2 * n); i++) {
		DiagonalOmegaC[i][i] = OmegaC[i];
	}

	double ProccessNoise[6][6];
	Matrix_Zero(n, n, ProccessNoise);
	ProccessNoise[0][0] = 1e-35;
	ProccessNoise[1][1] = 1e-35;
	ProccessNoise[2][2] = 1e-35;
	ProccessNoise[3][3] = 3e-10;
	ProccessNoise[4][4] = 3e-10;
	ProccessNoise[5][5] = 3e-10;
	MeasurmentNoise = pow((100e-9), 2); //(T)^2   (1e-4 * (10e-6))
	ZHat = 0;
	Matrix_Zero(n, (2 * n) + 1, YY);
	VectorZeros((2 * n) + 1, GammaK);
	Matrix_Zero(n, (2 * n) + 1, Chi);
	Matrix_Zero(n, (2 * n) + 1, xHatMinus);
	Matrix_Zero(n, (2 * n) + 1, Y1);
	mu = 398600e+9; // m^3 * s^-2
	VectorZeros(3, EstimateOfR);
	VectorZeros(3, EstimateOfV);
}
void Predict()
{
	CalcSigma(xHat, PMinus, gamma, n, ChiMinus);
	OrbitDynamic();
	double xHatpp[6];
	VectorZeros(n, xHatpp);
	for (int iS = 0; iS < (2 * n + 1); iS++)
	{
		double x[6], xHatNew[6];
		VectorZeros(n, xHatNew);
		VectorZeros(n, x);
		double SigmaFromRealDynamics[6];

		for (int j = 0; j < n; j++)
		{
			x[j] = ChiMinus[j][iS];
		}
		RungeKutta4(0, n, x, 0.25, 1, SigmaFromRealDynamics);
		for (int j = 0; j < n; j++)
		{
			Chi[j][iS] = SigmaFromRealDynamics[j];
		}

		Vector_Multiplier2Num(n, SigmaFromRealDynamics, OmegaS[iS], xHatNew);
		Vector_Sum(n, xHatpp, xHatNew, xHatpp);

	}
	for (int i = 0; i < (2 * n + 1); i++)
	{
		for (int j = 0; j < n; j++)
		{
			xHatMinus[j][i] = xHatpp[j];
		}
	}

	Matrix_Subtract6to13(n, (2 * n + 1), Chi, xHatMinus, Y1);


	double Y1Trans[13][6];
	Matrix_Transpose6to13(n, (2 * n + 1), Y1, Y1Trans);

	double diagocy1[13][6];
	double diagocy2[6][6];
	Matrix_MultiplierAsymmetric13to6((2 * n + 1), (2 * n + 1), DiagonalOmegaC, (2 * n + 1), n, Y1Trans, diagocy1);
	Matrix_MultiplierAsymmetric6to6(n, (2 * n + 1), Y1, (2 * n + 1), n, diagocy1, diagocy2);

	Matrix_Sum6to6(n, n, diagocy2, ProccessNoise, PMinusP);

	EqualVector(6, xHatpp, xHatP);

}
void Update(double NormMagnetInNED[], double InputCurrentMeasurment)
{
	CalcSigma(xHatP, PMinusP, gamma, n, ChiPlus);
	ZHat = ZHat * 0;
	for (int iS = 0; iS < (2 * n) + 1; iS++)
	{
		GammaK[iS] = NormMagnetInNED[iS];
		ZHat += OmegaS[iS] * NormMagnetInNED[iS];
	}

	double ZHatA[(2 * 6) + 1];
	for (int i = 0; i < (2 * n) + 1; i++) {
		ZHatA[i] = ZHat;
	}

	Vector_Subtract((2 * n) + 1, GammaK, ZHatA, Y2);
	double MutY2Diag[(2 * 6) + 1];
	Matrix13to13_Multiplier2Vector13((2 * n) + 1, DiagonalOmegaC, Y2, MutY2Diag);
	double DotY = 0;
	DotY = Vector_Dot((2 * n) + 1, Y2, MutY2Diag);
	Pzz = DotY + MeasurmentNoise;

	//Matrix_Subtract6to13(n, (2 * n) + 1, ChiPlus, xHatMinus, Y3);
	Matrix6to13_Multiplier2Vector13(n, Y1, MutY2Diag, Pxz);

	Vector_Divide2Num(n, Pxz, Pzz, k);
	double KK[6];
	Vector_Multiplier2Num(n, k, (InputCurrentMeasurment - ZHat), KK);
	Vector_Sum(n, xHatP, KK, xHat);

	double kMatrixV[6][1];
	double kMatrixH[1][6];
	double KVH[6][6];
	for (int i = 0; i < n; i++) {
		kMatrixV[i][0] = k[i] * Pzz;
		kMatrixH[0][i] = k[i];
	}
	Matrix_MultiplierAsymmetric6116(n, 1, kMatrixV, 1, n, kMatrixH, KVH);
	Matrix_Subtract6to6(n, n, PMinusP, KVH, PMinus);

}
void CalcSigma(double state[], double covariance[6][6], double gamma, int n, double Chi_Output[6][13])
{
	/*Chi_Output = (double**)malloc(sizeof(double*) * n);
	for (int i = 0; i < n; i++)
	{
		Chi_Output[i] = (double*)malloc(sizeof(double) * (2 * n + 1));
	}*/
	double Chol[6][6];
	Matrix_Zero6(Chol);
	cholesky6to6(6, covariance, Chol);
	//double* CholVec;// = (double*)calloc(n * n, sizeof(double));
	////double* CovVec = (double*)calloc(n * n, sizeof(double));
	//double CovVec[36];
	//int kkk = 0;
	//for (int i = 0; i < n; i++) 
	//{
	//	for (int j = 0; j < n; j++)
	//	{
	//		CovVec[kkk] = covariance[i][j];
	//		kkk++;
	//	}
	//}
	//
	////cholesky6to6(n, covariance, Chol);
	//cholesky(&CovVec, n, CholVec);
	//for (int i = 0; i < n; i++)
	//{
	//	for (int j = 0; j < n; j++)
	//	{
	//		Chol[i][j] = *(CholVec + i + j);
	//	}
	//}
	Matrix_Multiplier2Num6to6(n, n, Chol, gamma, A);
	double Y[6][6];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Y[i][j] = state[i];
		}
	}
	double YplusA[6][6];
	double YminusA[6][6];
	Matrix_Sum6to6(n, n, Y, A, YplusA);
	Matrix_Subtract6to6(n, n, Y, A, YminusA);

	for (int i = 0; i < n; i++)
	{
		Chi_Output[i][0] = state[i];
	}
	for (int j = 1; j < n + 1; j++)
	{
		for (int i = 0; i < n ; i++)
		{
			Chi_Output[i][j] = YplusA[i][j-1];
		}
	}
	for (int j = n + 1; j < (2 * n) + 1; j++)
	{
		for (int i = 0; i < n; i++)
		{
			Chi_Output[i][j] = YminusA[i][j-(n + 1)];
		}
	}
}
void RungeKutta4(double t0, int lengthOfState, double State[], double h, double t1, double Output[])
{
	double DeltaT = t1 - t0;
	int n = (int)(DeltaT / h);
	double k1[6], k12[6], k1p2[6], k1p21[6], kAll[6];
	double k2[6], k2t2[6], k2p2[6];
	double k3[6], k3t2[6];
	double k4[6], k34[6];

	double y[6], y1[6], y2[6], y3[6];
	double initial[6];
	for (int j = 0; j < lengthOfState; j++)
	{
		y[j] = State[j];
		initial[j] = State[j];

	}


		// Step1
		dynamic(t0, y, InputQuatOfSatECIToBody, k1);
		//Vector_Multiplier2Num(lengthOfState, k1, h, k1);
		// Step2
		Vector_Multiplier2Num(lengthOfState, k1, 0.5, k12);
		Vector_Sum(lengthOfState, initial, k12, y1);
		dynamic((t0 + 0.5 * h), y1, InputQuatOfSatECIToBody, k2);
		// Step3
		Vector_Multiplier2Num(lengthOfState, k2, 0.5, k2t2);
		Vector_Sum(lengthOfState, initial, k2t2, y2);
		dynamic((t0 + 0.5 * h), y2, InputQuatOfSatECIToBody, k3);
		// Step4
		Vector_Sum(lengthOfState, initial, k3, y3);
		dynamic((t0 + h), y3, InputQuatOfSatECIToBody, k4);
		// end stepping
		Vector_Multiplier2Num(lengthOfState, k2, 2, k2p2);
		Vector_Sum(lengthOfState, k1, k2p2, k1p2);
		Vector_Multiplier2Num(lengthOfState, k3, 2, k3t2);
		Vector_Sum(lengthOfState, k3t2, k4, k34);
		Vector_Sum(lengthOfState, k34, k1p2, k1p21);
		Vector_Multiplier2Num(lengthOfState, k1p21, (DeltaT / 6.0), kAll);
		Vector_Sum(lengthOfState, initial, kAll, Output);

		t0 = t0 + h;

}

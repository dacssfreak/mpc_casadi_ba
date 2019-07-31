#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"
#include "main.h"

//define dimensions of the system
int rowA = 4, colA = 4;
int rowx = 4; colx = 1;
int rowu = 1; colu = 1;
int rowC = 1; colC = 4;
int rowR = 1; colR = 1;

//define needed matrices
Matrix* x_estd;
Matrix* u;
Matrix* P;
Matrix* Q;
Matrix* C;
Matrix* R;
Matrix* I;
Matrix* x_trued;
Matrix* V;

FILE* fp;
FILE* gnuX;
FILE* gnuY;
FILE* gnuZ;
FILE* gnuQ;

int samples = 1000;
float delta = 0.01f;		//sample time
float startTime = 0.0f;

float currentTime;
float currentPosition;
float currentVelocity;
float currentAngle;
float currentAngularVelocity;

float currentPositionEst;
float currentVelocityEst;
float currentAngleEst;
float currentAngularVelocityEst;

float randn(float mu, float asd)
{
	float U1, U2, W, mult;
	float X1;

	do
	{
		U1 = -1.0f + ((float)rand() / RAND_MAX) * 2.0f;
		U2 = -1.0f + ((float)rand() / RAND_MAX) * 2.0f;
		W = (float)pow(U1, 2) + (float)pow(U2, 2);
	} while (W >= 1 || W == 0);

	mult = sqrt((-2.0f * log(W)) / W);
	X1 = U1 * mult;

	return mu + asd * X1;
}

Matrix* KalmanCalculations(int rep)
{
	//update A (u = 0)
	Matrix* A = Calc_A();

	//discretize A with Tustin transform
	Matrix* A_temp = Calc_A_temp(A);
	
	// A_d = exp(A*delta) = (I + 0.5 * A * delta) * (I - 0.5 * A * delta)^-1;	A_temp = 0.5 * A * delta
	Matrix* IA2inv = matrix_alloc(rowA, colA);
	Matrix* IsubA_temp = matrix_subtract(I, A_temp);
	matrix_inverse(IsubA_temp, IA2inv);
	matrix_free(IsubA_temp);
	Matrix* IaddA_temp = matrix_add(I, A_temp);
	Matrix* A_d = matrix_multiply(IaddA_temp, IA2inv);
	matrix_free(IaddA_temp);
	matrix_free(IA2inv);
		
	//update W
	Matrix* W = Calc_W();
	
	//true state
	Matrix* x_truec = Calc_x_truec();

	//calc x trued
	x_trued = RK4(delta, x_truec, x_trued);

	matrix_print_detail(x_trued, "x_trued", rep+1);

	//measurement noise
	Matrix* v = matrix_alloc(rowC, colx);
	v->matrix_entry[0][0] = randn(0.0f, 0.01f);

	//project the state ahead
	Matrix* x_estc = Calc_x_estc();

	x_estd = RK4(delta, x_estc, x_estd);

	//Project the error covariance ahead, P = A * P * A' + W * Q * W'
	Matrix* A_dP = matrix_multiply(A_d, P);
	Matrix* WQ = matrix_multiply(W, Q);
	Matrix* A_dt = matrix_transpose(A_d);
	Matrix* Wt = matrix_transpose(W);
	Matrix* A_dPA_dt = matrix_multiply(A_dP, A_dt);
	Matrix* WQWt = matrix_multiply(WQ, Wt);
	matrix_free(P);
	P = matrix_add(A_dPA_dt, WQWt);
	matrix_free(A_dP);
	matrix_free(WQ);
	matrix_free(A_dt);
	matrix_free(Wt);
	matrix_free(A_dPA_dt);
	matrix_free(WQWt);

	//K = P * C' * (C * P * C' + V * R * V')^-1; compute the kalman gain
	Matrix* CPCt_VRVt_inv = matrix_alloc(rowC, rowC);
	//C*P*C' + V*R*V'
	Matrix* Ct = matrix_transpose(C);
	Matrix* PCt = matrix_multiply(P, Ct);
	Matrix* CPCt = matrix_multiply(C, PCt);
	Matrix* Vt = matrix_transpose(V);
	Matrix* VR = matrix_multiply(V, R);
	Matrix* VRVt = matrix_multiply(VR, Vt);
	Matrix* CPCt_VRVt = matrix_add(CPCt, VRVt);
	matrix_free(CPCt);
	matrix_free(Vt);
	matrix_free(VR);
	matrix_free(VRVt);
	//matrix_inverse(CPCtR, CPCtRinv);			//matrix_inverse only works for 2 or higher dimensions
	CPCt_VRVt_inv->matrix_entry[0][0] = 1 / CPCt_VRVt->matrix_entry[0][0];
	Matrix* K_kf = matrix_multiply(PCt, CPCt_VRVt_inv);
	matrix_free(Ct);
	matrix_free(PCt);
	matrix_free(CPCt_VRVt);
	matrix_free(CPCt_VRVt_inv);

	//update estimate with measurement, x_est = x_est + K * (y_mes - y_est)
	Matrix* Cx_trued = matrix_multiply(C, x_trued);
	Matrix* Cx_estd = matrix_multiply(C, x_estd);
	Matrix* Cx_truedaddv = matrix_add(Cx_trued, v);
	Matrix* delta_y = matrix_subtract(Cx_truedaddv, Cx_estd);
	Matrix* K_kfdelta_y = matrix_multiply(K_kf, delta_y);
	Matrix* x_estd_temp = matrix_add(x_estd, K_kfdelta_y);
	matrix_free(x_estd);
	x_estd = x_estd_temp;
	matrix_free(Cx_trued);
	matrix_free(Cx_estd);
	matrix_free(Cx_truedaddv);
	matrix_free(delta_y);
	matrix_free(K_kfdelta_y);
	matrix_print_detail(x_estd, "x_estd", rep + 1);

	//update error covariance, P = (I - K * C) * P
	Matrix* K_kfC = matrix_multiply(K_kf, C);
	Matrix* IsubK_kfC = matrix_subtract(I, K_kfC);
	Matrix* P_temp2 = matrix_multiply(IsubK_kfC, P);
	matrix_free(K_kfC);
	matrix_free(IsubK_kfC);
	matrix_free(P);
	P = P_temp2;

	//check convergence, x_est - x_true
	Matrix* delta_x = matrix_subtract(x_estd, x_trued);
	matrix_print_detail(delta_x, "delta_x", rep+1);

	//GnuPlotStuff1(rep);

	matrix_free(A);
	matrix_free(A_temp);
	matrix_free(A_d);
	matrix_free(W);
	matrix_free(x_truec);
	matrix_free(v);
	matrix_free(x_estc);
	matrix_free(K_kf);
	matrix_free(delta_x);
	
	return x_estd;
}

void GnuPlotStuff1(int rep)
{
	currentPosition = x_trued->matrix_entry[0][0];
	currentVelocity = x_trued->matrix_entry[1][0];
	currentAngle = x_trued->matrix_entry[2][0];
	currentAngularVelocity = x_trued->matrix_entry[3][0];

	currentPositionEst = x_estd->matrix_entry[0][0];
	currentVelocityEst = x_estd->matrix_entry[1][0];
	currentAngleEst = x_estd->matrix_entry[2][0];
	currentAngularVelocityEst = x_estd->matrix_entry[3][0];

	/* open the file for writing*/
	fp = fopen("C:\\Users\\Faker\\source\\repos\\EKFc\\EKFc\\test.txt", "a");

	fprintf(fp, "\n%f \t %f \t %f \t %f \t\t %f \t %f \t %f \t %f \t %f", currentTime, currentPosition, currentVelocity, currentAngle, currentAngularVelocity, currentPositionEst, currentVelocityEst, currentAngleEst, currentAngularVelocityEst);

	/* close the file*/
	fclose(fp);

	currentTime = currentTime + delta;

	if (rep % 30 == 0) {
		fprintf(gnuX, "plot 'C:\\Users\\Faker\\source\\repos\\EKFc\\EKFc\\test.txt' using 1:2 title 'True Position', 'C:\\Users\\Faker\\source\\repos\\EKFc\\EKFc\\test.txt' using 1:6 title 'Estimated Position' \n");
		fprintf(gnuY, "plot 'C:\\Users\\Faker\\source\\repos\\EKFc\\EKFc\\test.txt' using 1:3 title 'True Velocity', 'C:\\Users\\Faker\\source\\repos\\EKFc\\EKFc\\test.txt' using 1:7 title 'Estimated Velocity' \n");
		fprintf(gnuZ, "plot 'C:\\Users\\Faker\\source\\repos\\EKFc\\EKFc\\test.txt' using 1:4 title 'True Angle', 'C:\\Users\\Faker\\source\\repos\\EKFc\\EKFc\\test.txt' using 1:8 title 'Estimated Angle' \n");
		fprintf(gnuQ, "plot 'C:\\Users\\Faker\\source\\repos\\EKFc\\EKFc\\test.txt' using 1:5 title 'True Angular Velocity', 'C:\\Users\\Faker\\source\\repos\\EKFc\\EKFc\\test.txt' using 1:9 title 'Estimated Angular Velocity' \n");

		fflush(gnuY);
		fflush(gnuX);
		fflush(gnuZ);
		fflush(gnuQ);
	}
}

Matrix* Calc_W()
{
	Matrix* W = matrix_alloc(rowA, colA);

	W->matrix_entry[0][0] = x_estd->matrix_entry[0][0];
	W->matrix_entry[0][1] = 0.0f;
	W->matrix_entry[0][2] = 0.0f;
	W->matrix_entry[0][3] = 0.0f;
	W->matrix_entry[1][0] = 0.0f;
	W->matrix_entry[1][1] = x_estd->matrix_entry[1][0];
	W->matrix_entry[1][2] = 0.0f;
	W->matrix_entry[1][3] = 0.0f;
	W->matrix_entry[2][0] = 0.0f;
	W->matrix_entry[2][1] = 0.0f;
	W->matrix_entry[2][2] = x_estd->matrix_entry[2][0];
	W->matrix_entry[2][3] = 0.0f;
	W->matrix_entry[3][0] = 0.0f;
	W->matrix_entry[3][1] = 0.0f;
	W->matrix_entry[3][2] = 0.0f;
	W->matrix_entry[3][3] = x_estd->matrix_entry[3][0];

	return W;
}

Matrix* Calc_A_temp(Matrix* A)
{
	Matrix* A_temp = matrix_alloc(rowA, colA);

	A_temp->matrix_entry[0][0] = (0.5f * A->matrix_entry[0][0] * delta);
	A_temp->matrix_entry[0][1] = (0.5f * A->matrix_entry[0][1] * delta);
	A_temp->matrix_entry[0][2] = (0.5f * A->matrix_entry[0][2] * delta);
	A_temp->matrix_entry[0][3] = (0.5f * A->matrix_entry[0][3] * delta);
	A_temp->matrix_entry[1][0] = (0.5f * A->matrix_entry[1][0] * delta);
	A_temp->matrix_entry[1][1] = (0.5f * A->matrix_entry[1][1] * delta);
	A_temp->matrix_entry[1][2] = (0.5f * A->matrix_entry[1][2] * delta);
	A_temp->matrix_entry[1][3] = (0.5f * A->matrix_entry[1][3] * delta);
	A_temp->matrix_entry[2][0] = (0.5f * A->matrix_entry[2][0] * delta);
	A_temp->matrix_entry[2][1] = (0.5f * A->matrix_entry[2][1] * delta);
	A_temp->matrix_entry[2][2] = (0.5f * A->matrix_entry[2][2] * delta);
	A_temp->matrix_entry[2][3] = (0.5f * A->matrix_entry[2][3] * delta);
	A_temp->matrix_entry[3][0] = (0.5f * A->matrix_entry[3][0] * delta);
	A_temp->matrix_entry[3][1] = (0.5f * A->matrix_entry[3][1] * delta);
	A_temp->matrix_entry[3][2] = (0.5f * A->matrix_entry[3][2] * delta);
	A_temp->matrix_entry[3][3] = (0.5f * A->matrix_entry[3][3] * delta);

	return A_temp;
}

Matrix* Calc_A()
{
	Matrix* A = matrix_alloc(rowA, colA);
	A->matrix_entry[0][0] = 0.0f;
	A->matrix_entry[0][1] = 1.0f;
	A->matrix_entry[0][2] = 0.0f;
	A->matrix_entry[0][3] = 0.0f;
	A->matrix_entry[1][0] = 0.0f;
	A->matrix_entry[1][1] = 0.0f;
	A->matrix_entry[1][2] = (float)(((2.0f * pow(x_estd->matrix_entry[3][0], 2) * cos(x_estd->matrix_entry[2][0]) + 20.0f * pow(cos(x_estd->matrix_entry[2][0]), 2) - 10.0f) / (pow(cos(x_estd->matrix_entry[2][0]), 2) - 6.0f)) + ((sin(2.0f * x_estd->matrix_entry[2][0]) * (2.0f * sin(x_estd->matrix_entry[2][0]) * pow(x_estd->matrix_entry[3][0], 2) + 5.0f * sin(2.0f * x_estd->matrix_entry[2][0]))) / (pow((pow(cos(x_estd->matrix_entry[2][0]), 2) - 6.0f), 2))));
	A->matrix_entry[1][3] = (float)((4.0f * x_estd->matrix_entry[3][0] * sin(x_estd->matrix_entry[2][0])) / (pow(cos(x_estd->matrix_entry[2][0]), 2) - 6));
	A->matrix_entry[2][0] = 0.0f;
	A->matrix_entry[2][1] = 0.0f;
	A->matrix_entry[2][2] = 0.0f;
	A->matrix_entry[2][3] = 1.0f;
	A->matrix_entry[3][0] = 0.0f;
	A->matrix_entry[3][1] = 0.0f;
	A->matrix_entry[3][2] = (float)(((60.0f * cos(x_estd->matrix_entry[2][0]) + 2.0f * pow(x_estd->matrix_entry[3][0], 2) * (2.0f * pow(cos(x_estd->matrix_entry[2][0]), 2) - 1.0f)) / (2.0f * pow(cos(x_estd->matrix_entry[2][0]), 2) - 12.0f)) + ((2.0f * sin(2.0f * x_estd->matrix_entry[2][0]) * (sin(2.0f * x_estd->matrix_entry[2][0]) * pow(x_estd->matrix_entry[3][0], 2) + 60.0f * sin(x_estd->matrix_entry[2][0]))) / (pow((2.0f * pow(cos(x_estd->matrix_entry[2][0]), 2) - 12.0f), 2))));
	A->matrix_entry[3][3] = (float)((2.0f * x_estd->matrix_entry[3][0] * sin(2.0f * x_estd->matrix_entry[2][0])) / (2.0f * pow(cos(x_estd->matrix_entry[2][0]), 2) - 12.0f));

	return A;
}

Matrix* Calc_x_estc()
{
	Matrix* x_estc = matrix_alloc(rowx, colx);

	x_estc->matrix_entry[0][0] = (float)x_estd->matrix_entry[1][0];
	x_estc->matrix_entry[1][0] = (float)((u->matrix_entry[0][0] - 2.0f * sin(x_estd->matrix_entry[2][0]) * pow(x_estd->matrix_entry[3][0], 2) - 5.0f * sin(2.0f * x_estd->matrix_entry[2][0])) / (6.0f - pow(cos(x_estd->matrix_entry[2][0]), 2)));
	x_estc->matrix_entry[2][0] = (float)x_estd->matrix_entry[3][0];
	x_estc->matrix_entry[3][0] = (float)((u->matrix_entry[0][0] * cos(x_estd->matrix_entry[2][0]) - 60.0f * sin(x_estd->matrix_entry[2][0]) - sin(2.0f * x_estd->matrix_entry[2][0]) * pow(x_estd->matrix_entry[3][0], 2)) / (12.0f - 2.0f * pow(cos(x_estd->matrix_entry[2][0]), 2)));

	return x_estc;
}

Matrix* Calc_x_truec()
{
	Matrix* x_truec = matrix_alloc(rowx, colx);

	x_truec->matrix_entry[0][0] = (float)(x_trued->matrix_entry[1][0] + randn(0.0f, 0.01f) * x_trued->matrix_entry[0][0]);
	x_truec->matrix_entry[1][0] = (float)(((u->matrix_entry[0][0] - 2.0f * sin(x_trued->matrix_entry[2][0]) * pow(x_trued->matrix_entry[3][0], 2) - 5.0f * sin(2.0f * x_trued->matrix_entry[2][0])) / (6.0f - pow(cos(x_trued->matrix_entry[2][0]), 2))) + randn(0.0f, 0.01f) * x_trued->matrix_entry[1][0]);
	x_truec->matrix_entry[2][0] = (float)(x_trued->matrix_entry[3][0] + randn(0.0f, 0.01f) * x_trued->matrix_entry[2][0]);
	x_truec->matrix_entry[3][0] = (float)(((u->matrix_entry[0][0] * cos(x_trued->matrix_entry[2][0]) - 60.0f * sin(x_trued->matrix_entry[2][0]) - sin(2.0f * x_trued->matrix_entry[2][0]) * pow(x_trued->matrix_entry[3][0], 2)) / (12.0f - 2.0f * pow(cos(x_trued->matrix_entry[2][0]), 2))) + randn(0.0f, 0.01f) * x_trued->matrix_entry[3][0]);
	
	return x_truec;
}

int main()
{
	//initial condition for x_true
	x_trued = matrix_alloc(rowx, colx);
	Fill_x_trued();
	
	//initial estimation
	x_estd = matrix_alloc(rowx, colx);
	Fill_x_estd();

	//input
	u = matrix_alloc(rowu, colu);
	Fill_u();

	//output matrix
	C = matrix_alloc(rowC, colC);
	Fill_C();

	//process noise matrix
	Q = matrix_alloc(rowA, colA);
	Fill_Q();

	//measurement error covariance matrix
	R = matrix_alloc(rowR, colR);
	R->matrix_entry[0][0] = 0.01f;

	//error covariance
	P = matrix_alloc(rowA, colA);
	Fill_P();

	//jacobian matrix of partial derivatives of measurement function with respect to v
	V = matrix_alloc(rowR, colR);
	Fill_V();
	
	I = matrix_callalloc(rowA);

	//GnuPlot_InitValues();

	for (int j = 0; j < samples; j++)
	{
		KalmanCalculations(j);
		print_alloc();
	}

	getch();
	
	return 0;
}

void GnuPlot_InitValues()
{
	gnuX = _popen("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persist", "w");
	fprintf(gnuX, "set xlabel 'Time [s]'\n");
	fprintf(gnuX, "set ylabel 'Position [m]'\n");

	gnuY = _popen("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persist", "w");
	fprintf(gnuY, "set xlabel 'Time [s]'\n");
	fprintf(gnuY, "set ylabel 'Velocity [m/s]'\n");

	gnuZ = _popen("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persist", "w");
	fprintf(gnuZ, "set xlabel 'Time [s]'\n");
	fprintf(gnuZ, "set ylabel 'Angle [rad]'\n");

	gnuQ = _popen("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persist", "w");
	fprintf(gnuQ, "set xlabel 'Time [s]'\n");
	fprintf(gnuQ, "set ylabel 'Angular Velocity [rad/s]'\n");

	float startPosition = x_trued->matrix_entry[0][0];
	float startVelocity = x_trued->matrix_entry[1][0];
	float startAngle = x_trued->matrix_entry[2][0];
	float startAngularVelocity = x_trued->matrix_entry[3][0];

	float startPositionEst = x_estd->matrix_entry[0][0];
	float startVelocityEst = x_estd->matrix_entry[1][0];
	float startAngleEst = x_estd->matrix_entry[2][0];
	float startAngularVelocityEst = x_estd->matrix_entry[3][0];

	/* open the file for writing*/
	fp = fopen("C:\\Users\\Faker\\source\\repos\\EKFc\\EKFc\\test.txt", "w");

	fprintf(fp, "Time \t\t Position \t Velocity \t Angle \t\t Angular Velocity \t PositionEst \t VelocityEst \t AngleEst \t Angular VelocityEst\n");
	fprintf(fp, "%f \t %f \t %f \t %f \t\t %f \t %f \t %f \t %f \t %f", startTime, startPosition, startVelocity, startAngle, startAngularVelocity, startPositionEst, startVelocityEst, startAngleEst, startAngularVelocityEst);

	/* close the file*/
	fclose(fp);
	currentTime = startTime + delta;
}

void Fill_V()
{
	V->matrix_entry[0][0] = 1.0f;
}

void Fill_P()
{
	P->matrix_entry[0][0] = 10.0f;
	P->matrix_entry[0][1] = 0.0f;
	P->matrix_entry[0][2] = 0.0f;
	P->matrix_entry[0][3] = 0.0f;
	P->matrix_entry[1][0] = 0.0f;
	P->matrix_entry[1][1] = 10.0f;
	P->matrix_entry[1][2] = 0.0f;
	P->matrix_entry[1][3] = 0.0f;
	P->matrix_entry[2][0] = 0.0f;
	P->matrix_entry[2][1] = 0.0f;
	P->matrix_entry[2][2] = 10.0f;
	P->matrix_entry[2][3] = 0.0f;
	P->matrix_entry[3][0] = 0.0f;
	P->matrix_entry[3][1] = 0.0f;
	P->matrix_entry[3][2] = 0.0f;
	P->matrix_entry[3][3] = 10.0f;
}

void Fill_Q()
{
	Q->matrix_entry[0][0] = 0.1f;
	Q->matrix_entry[0][1] = 0.0f;
	Q->matrix_entry[0][2] = 0.0f;
	Q->matrix_entry[0][3] = 0.0f;
	Q->matrix_entry[1][0] = 0.0f;
	Q->matrix_entry[1][1] = 0.1f;
	Q->matrix_entry[1][2] = 0.0f;
	Q->matrix_entry[1][3] = 0.0f;
	Q->matrix_entry[2][0] = 0.0f;
	Q->matrix_entry[2][1] = 0.0f;
	Q->matrix_entry[2][2] = 0.1f;
	Q->matrix_entry[2][3] = 0.0f;
	Q->matrix_entry[3][0] = 0.0f;
	Q->matrix_entry[3][1] = 0.0f;
	Q->matrix_entry[3][2] = 0.0f;
	Q->matrix_entry[3][3] = 0.1f;
}

void Fill_C()
{
	C->matrix_entry[0][0] = 0.0f;
	C->matrix_entry[0][1] = 1.0f;
	C->matrix_entry[0][2] = 0.0f;
	C->matrix_entry[0][3] = 0.0f;
}

void Fill_u()
{
	u->matrix_entry[0][0] = 0.0f;
}

void Fill_x_estd()
{
	x_estd->matrix_entry[0][0] = 0.0f;
	x_estd->matrix_entry[1][0] = 0.0f;
	x_estd->matrix_entry[2][0] = 1.3f;
	x_estd->matrix_entry[3][0] = 0.0f;
}

void Fill_x_trued()
{
	x_trued->matrix_entry[0][0] = 0.0f;
	x_trued->matrix_entry[1][0] = 0.0f;
	x_trued->matrix_entry[2][0] = 1.0f;
	x_trued->matrix_entry[3][0] = 0.0f;
}

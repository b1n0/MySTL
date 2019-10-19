#include "h.h"
#define ALPHA 0.
#define BETA 1.

double u(double* y);
double u(double* y) { return y[3] > 0 ? 2. : -2.; }

void f(double x, double* y, int size, double* res) {
	res[0] = u(y);
	res[1] = y[0];
	res[2] = -BETA*pow(y[1], BETA - 1)/(1 + exp(ALPHA*y[0]*y[0]));
	res[3] = -1*y[2] + 2*ALPHA*pow(y[1], BETA)*y[0]*exp(ALPHA*y[0]*y[0])/pow(1+exp(ALPHA*y[0]*y[0]), 2);
}

double eigen_value(double x, double* y) { return 0; }

void discrepancy(double* y0, double* y, double* v) { 
	v[0] = y[0];
	v[1] = y[1] + y0[1];
	v[2] = y[2] + y0[2];	
}

void jacobian(double** m, double* y0, double* v) {
	double y[4], dv[3];
	for(int i = 0; i < 3; i++) {
		y0[1+i] += DELTA;
		runge_hardcore(0, 4, y0, y, 4, 1.e-8, 1.e-7);
		y0[1+i] -= DELTA;
		discrepancy(y0, y, dv);
		m[0][i] = (dv[0] - v[0])/DELTA;
		m[1][i] = (dv[1] - v[1])/DELTA;
		m[2][i] = (dv[2] - v[2])/DELTA;
	}
}

int main(void) {
	double a, b, y0[4], y[4];
	y0[0] = 0;
	shoot(a, b, y0, 4, 1, 0.0001);
	printf("%lf %lf %lf %lf \n", y0[0], y0[1], y0[2], y0[3]);	
	return 0;
}

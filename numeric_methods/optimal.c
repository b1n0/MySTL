#include "h.h"
#define ALPHA 0.
#define BETA 1.

double u(double* y);
double u(double* y) { return y[3] > 0 ? 2. : -2.; }

void start_value(double* y0) { y0[1] = 0.; y0[2] = 1.0; y0[3] = 3./4.; }

double eigen_value(double x, double* y) { x++; y++; return 0; }

void f(double x, double* y, double* res) {
	x++;
	res[0] = u(y); // y
	res[1] = y[0]; // x
	res[2] = -BETA*pow(y[1], BETA - 1)/(1 + exp(ALPHA*y[0]*y[0])); //px
	res[3] = -y[2] + 2*ALPHA*pow(y[1], BETA)*y[0]*exp(ALPHA*y[0]*y[0])/pow(1+exp(ALPHA*y[0]*y[0]), 2); //py
}
void discrepancy(double* y0, double* y, double* res) { 
	res[0] = y[0];
	res[1] = y[1] + y0[1];
	res[2] = y[2] + y0[2];	
}

double horde_method(double x1, double x2, double* y1, double* y2, int size) { 
	double y[ST_SIZE], xroot;
	while(y1[3]*y2[3] < 0. && (fabs(y1[3]) > 1.e-10 || fabs(y2[3]) > 1.e-10)) {
		xroot = x1 - y1[3]*(x2 - x1)/(y2[3] - y1[3]);
		dormand8(x1, y1, y, size, xroot - x1);
		if(y1[3]*y[3] > 0) { memcpy(y1, y, size*sizeof(double)); x1 = xroot; }
		else { memcpy(y2, y, size*sizeof(double)); x2 = xroot; }
	}
	return x2;
}

int main(void) {
	double a, b, y0[4];
	a = 0.; b = 4.;
	y0[0] = 0;
	if(shoot(a, b, y0, 4, 1, 0.001, discrepancy) == 0) {
		printf("%lf %lf %lf %lf \n", y0[0], y0[1], y0[2], y0[3]);	
		track(a, b, y0, 4, 100, 1);
	}
	return 0;
}

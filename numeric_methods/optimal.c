#include "h.h"
#define ALPHA 0.
#define BETA 1.

double u(double* y) { return y[4] > 0 ? 2. : -2.; }
double g(double* y) { return y[4]; }
void start_value(double* y0)  { y0[2] = 0.; y0[3] = 1.0; y0[4] = 3./4.; }
double eigen_value(double x, double* y) { x++; y++; return 1.; }

void f(double x, double* y, double* res) {
	x++;
	res[0] = pow(y[2], BETA)/(1. + exp(ALPHA*y[1]*y[1]));
	res[1] = u(y); // y
	res[2] = y[1]; // x
	res[3] = -BETA*pow(y[2], BETA - 1)/(1 + exp(ALPHA*y[1]*y[1])); //px
	res[4] = -y[3] + 2*ALPHA*pow(y[2], BETA)*y[1]*exp(ALPHA*y[1]*y[1])/pow(1+exp(ALPHA*y[1]*y[1]), 2); //py
}
void discrepancy(double* y0, double* y, double* res) { 
	res[0] = y[1];
	res[1] = y[2] + y0[2];
	res[2] = y[3] + y0[3];	
}

int main(void) {
	double a, b, y0[5];
	a = 0.; b = 4.;
	y0[0] = y0[1] = 0;
	if(shoot(a, b, y0, 5, 2, 0.001, discrepancy, 1.e-13, 1.e-12) == 0) {
		printf("%lf %lf %lf %lf %lf\n", y0[0], y0[1], y0[2], y0[3], y0[4]);	
		track(a, b, y0, 5, 100, 1);
		runge_numbers(a, b, y0, 5);
		printf("global error = %1.30lf \n", integrate_autostep(a, b, y0, y0, 5, 1.e-13, 1.e-12, 1.e-2, dormand8, 0));
		printf("B0 = %lf \n", y0[0]);
	}
	return 0;
}

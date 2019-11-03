#include "h.h"
#define NUM_POINTS 100
#define ALPHA 1.

void f(double x, double* y, double* res) {
	res[1] = (y[0] - y[1])/ALPHA;
	res[0] = (y[0] - y[1])/ALPHA - ALPHA*x/y[1];
}
double eigen_value(double x, double* y) { return sqrt(4/ALPHA - ALPHA*ALPHA*x*x/(y[1]*y[1])); }
void start_value(double* y0) { y0[0] = 0.; y0[1] = 2.; }
void discrepancy(double* y0, double* y, double* v) { v[0] = y[1]; }

int main(void) {
	double y0[2], a = 0., b = 1.;
	y0[0] = 0.; y0[1] = 0.;
	if(shoot(b, a, y0, 2, 1, 1.e-6, discrepancy) == 0) {
		printf("%lf %lf\n", y0[0], y0[1]);
		printf("global error = %5.20lf \n", track(b, a, y0, 2, NUM_POINTS));	
		runge_numbers(b, 0.75, y0, 2);
		runge_numbers(b, 0.5, y0, 2);
		runge_numbers(b, 0.25, y0, 2);
		runge_numbers(b, 0., y0, 2);
	}
	return 0;	
}

#include "h.h"
#define NUM_POINTS 100

void f(double x, double* y, double* res) {
	x++; 
	res[0] = y[1];
	res[1] = -y[0];
}
double eigen_value(double x, double* y) { x++; y++; return 0; }
double horde_method(double x1, double x2, double* y1, double* y2, int size) { return x2; }

int main(void) {
	int i = 0;
	double y0[2], y[2], err;
	y0[0] = 0.; y0[1] = 1.;
	for(double T = PI; i < 7; i++, T *= 10.) {
		err = integrate_autostep(0., T, y0, y, 2, 1.e-9, 1.e-8, 1.e-3, dormand8, 1);
		printf("%5.30lf | %5.30lf | %5.30lf | %lf |\n", fabs(y[0]), fabs(y[1] - cos(T)), err, T);
		err = integrate_autostep(0., T, y0, y, 2, 1.e-11, 1.e-10, 1.e-3, dormand8, 1);
		printf("%5.30lf | %5.30lf | %5.30lf | %lf |\n", fabs(y[0]), fabs(y[1] - cos(T)), err, T);
		err = integrate_autostep(0., T, y0, y, 2, 1.e-13, 1.e-12, 1.e-3, dormand8, 1);
		printf("%5.30lf | %5.30lf | %5.30lf | %lf |\n", fabs(y[0]), fabs(y[1] - cos(T)), err, T);
		runge_numbers(0., T, y0, 2);
		printf("\n");
	}	
	//printf("global error = %5.30lf\n", track(0, 2*PI, y0, 2, NUM_POINTS));
	return 0;
}

#include "h.h"
#define NUM_POINTS 100

void f(double x, double* y, double* res) {
	res[0] = y[1];
	res[1] = -y[0];
}
double eigen_value(double x, double* y) { return 0; }

int main(void) {
	double y[2], err;
	y[0] = 0.; y[1] = 1.;
	err = integrate_autostep(0., 2*PI, y, y, 2, 1.e-9, 1.e-8, 1.e-5);
	printf("global_error = %5.30lf\n%5.30lf\t%5.30lf\n", err, y[0], y[1]);
	y[0] = 0.; y[1] = 1.;
	runge_numbers(0., 2*PI, y, 2);
	printf("global error = %5.30lf\n", track(0, 2*PI, y, 2, NUM_POINTS));
	return 0;
}

#include "h.h"
#define NUM_POINTS 100

void f(double x, double* y, int size, double* res) {
	res[1] = y[0] - y[1];
	res[0] = y[0] - y[1] - x/y[1];
}

int main(void) {
	double y0[2];
	y0[0] = 0.; y0[1] = 0.;
        shooting(y0, 2, 1, 1., 0.);
	return 0;	
}

#include "h.h"
#define NUM_POINTS 100

void f(double x, double* y, int size, double* res) {
	res[0] = y[1] - y[0];
	res[1] = y[1] - y[0] - x/y[0];
}

int main(void) {
	double y0[5];
	y0[0] = 0.; y0[1] = 0;
        shooting(y0, 2, 1, 0., 1.);
	return 0;	
}

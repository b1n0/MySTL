#include "h.h"
#define NUM_POINTS 100

void f(double x, double* y, int size, double* res) {
	res[0] = y[1];
	res[1] = -1.*y[0];
}

int main(void) {
	double x0 = 0, x = 0, y0[2], y[2*NUM_POINTS], h = 2*3.14/NUM_POINTS;
	y0[0] = 0; y0[1] = 1;
	for(int i = 0; i < NUM_POINTS; i++, x+=h) {
		runge_with_autostep(x0, x, y0, y + 2*i, 2, 0.0001, 0.0000001, 10);
	}
	plot(y, NUM_POINTS);
	return 0;
}

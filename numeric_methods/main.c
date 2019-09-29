#include "h.h"
#define NUM_POINTS 100
#define delta 0.00000001

void f(double x, double* y, int size, double* res) {
	res[1] = y[0] - y[1];
	res[0] = y[0] - y[1] - x/y[1];
}

int main(void) {
	double y0[2], beta[2], y[2], m[2], h = 0.01, x = 1.;
	y0[0] = 0.; y0[1] = 0.;
	beta[0] = 0; beta[1] = 1.;
	while(1) {
		runge(1, 0, beta, y, 2, 10000);
		//runge_with_autostep(1, 0, beta, y, 2, 0.0001, 0.01);
		if (y[1]*y[1] < 0.00000001) break;
		printf("%lf\n", y[1]);
		beta[1] += delta;
		runge(1, 0, beta, m, 2, 10000);
		//runge_with_autostep(1, 0, beta, m, 2, 0.0001, 0.01);
		beta[1] -= delta;
		beta[1] += (-y[1])/((m[1]-y[1])/delta);
	}
	printf("%lf %lf \n", beta[1], y[1]);
	return 0;	
}


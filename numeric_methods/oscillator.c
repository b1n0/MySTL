#include "h.h"
#define NUM_POINTS 100

void f(double x, double* y, double* res) {
	res[0] = y[1];
	res[1] = -1.*y[0];
}
double eigen_value(double x, double* y) { return 0; }

int main(void) {
	double y0[2]; 
	y0[0] = 0; y0[1] = 1;

	printf("global error = %5.30lf\n", plot(0, 2*PI, y0, 2, NUM_POINTS));
		
	runge_numbers(0., 0.5*PI, y0, 2);
	runge_numbers(0., PI, y0, 2);
	runge_numbers(0., 1.5*PI, y0, 2);
	runge_numbers(0., 2*PI, y0, 2);
	
	return 0;
}

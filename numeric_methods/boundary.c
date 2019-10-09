#include "h.h"
#define NUM_POINTS 100
#define ALPHA 1.

void f(double x, double* y, int size, double* res) {
	size++;
	res[1] = (y[0] - y[1])/ALPHA;
	res[0] = (y[0] - y[1])/ALPHA - ALPHA*x/y[1];
}

double eigen_value(double x, double* y) { return sqrt(4/ALPHA - ALPHA*ALPHA*x*x/(y[1]*y[1])); }

int main(void) {
	double y0[2], y[2], u0[2], a = 0., b = 1., h = (b - a)/NUM_POINTS, err = 0.;
	FILE* f = fopen("track.txt", "w");
	y0[0] = 0.; y0[1] = 0.;
        
	shooting(b, a, y0, 2, 1, 1.e-6);
	printf("%lf %lf\n", y0[0], y0[1]);

	u0[0] = y0[0]; u0[1] = y0[1];
	for(double x = b + h, x0 = b; x > a; x -= h) {
		err = runge_with_autostep(x0, x-h, u0, y, 2, 1.e-9, 1.e-8);
		u0[0] = y[0]; u0[1] = y[1]; x0 = x - h;
		fprintf(f, "%lf %lf \n", y[0], y[1]);	
		//printf("%lf %lf\n", x, err);
	}
	fclose(f);
	printf("global error = %5.20lf \n", err);	
	plot("track.txt");
		
	runge_numbers(b, 0.75, y0, 2);
	/*
	runge_numbers(b, 0.5, y0, 2);
	runge_numbers(b, 0.25, y0, 2);
	runge_numbers(b, 0., y0, 2);
	*/
	return 0;	
}

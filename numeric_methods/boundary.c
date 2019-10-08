#include "h.h"
#define NUM_POINTS 100

void f(double x, double* y, int size, double* res) {
	size++;
	res[1] = (y[0] - y[1])/ALPHA;
	res[0] = (y[0] - y[1])/ALPHA - ALPHA*x/y[1];
}

int main(void) {
	double y0[2], y[2], a = 0., b = 1., h = (b - a)/NUM_POINTS, x = a, err = 0.;
	FILE* f = fopen("track.txt", "w");
	y0[0] = 0.; y0[1] = 0.;
        shooting(b, a, y0, 2, 1, 1.e-5);
	for(double x = b; x > a + h; x -= h) {
		err = runge(x, x-h, y0, y, 2, 1000);
		//err = runge_with_autostep(x, x-h, y0, y, 2, 1.e-8, 1.e-7);
		y0[0] = y[0]; y0[1] = y[1];
		fprintf(f, "%lf %lf \n", y[0], y[1]);	
	}
	printf("global error = %5.20lf \n", err);
	y0[0] = 0.; y0[1] = 0.;
	
	runge_numbers(b, 0.75, y0, 2);
	runge_numbers(b, 0.5, y0, 2);
	runge_numbers(b, 0.25, y0, 2);
	runge_numbers(b, 0., y0, 2);
	
	plot("track.txt");
	fclose(f);
	return 0;	
}

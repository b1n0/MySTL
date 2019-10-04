#include "h.h"
#define NUM_POINTS 10000
#define PI 3.14159265358979323846

void f(double x, double* y, int size, double* res) {
	size++; x++;
	res[0] = y[1];
	res[1] = -1.*y[0];
}

int main(void) {
	double x0 = 0, y0[2], y[2], h = 2*PI/NUM_POINTS, err = 0.;
	FILE *out = fopen("data.txt", "w");
	y0[0] = 0; y0[1] = 1;
	for(double x = x0; x < 2*PI - h; x+=h) {
		//euler(x0, x, y0, y, 2, 200000);
		//err = runge(x0, x, y0, y, 2, 2000);
		err = runge_with_autostep(x0, x, y0, y, 2, 1.e-9, 1.e-8);
		fprintf(out, "%lf %lf \n", y[0], y[1]);
	}
	printf("global error = %lf\n", err);
	runge_numbers(x0, 0.5*PI, y0, 2);
	runge_numbers(x0, PI, y0, 2);
	runge_numbers(x0, 1.5*PI, y0, 2);
	runge_numbers(x0, 2*PI, y0, 2);

	plot("data.txt");
	fclose(out);
	return 0;
}

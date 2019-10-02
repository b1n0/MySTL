#include "h.h"
#define NUM_POINTS 100
#define PI 3.14159265358979323846
void runge_numbers(double x0, double x, double* y0, int size);

void f(double x, double* y, int size, double* res) {
	res[0] = y[1];
	res[1] = -1.*y[0];
}

void plot(char* fname, int num_points) {
        FILE * gnuplot_pipe = popen("gnuplot -persistent", "w");
        fprintf(gnuplot_pipe, "%s '%s' %s\n", "plot", fname, "with line");
}

int main(void) {
	double x0 = 0, x = 0, y0[2], y[2], h = 2*PI/NUM_POINTS;
	y0[0] = 0; y0[1] = 1;
	FILE *out = fopen("data.txt", "w");
	for(int i = 0; i < NUM_POINTS; i++, x+=h) {
		runge_with_autostep(x0, x, y0, y, 2, 1.e-9, 1.e-8);
		fprintf(out, "%lf %lf \n", y[0], y[1]);
	}
	runge_numbers(x0, 0.5*PI, y0, 2);
	runge_numbers(x0, PI, y0, 2);
	runge_numbers(x0, 1.5*PI, y0, 2);
	runge_numbers(x0, 2*PI, y0, 2);

	plot("data.txt", NUM_POINTS);
	fclose(out);
	return 0;
}

void runge_numbers(double x0, double x, double* y0, int size) {
	double y8[ST_SIZE], y10[ST_SIZE], y12[ST_SIZE];
	runge_with_autostep(x0, x, y0, y8, 2, 1.e-9, 1.e-8);
	runge_with_autostep(x0, x, y0, y10, 2, 1.e-11, 1.e-10);
	runge_with_autostep(x0, x, y0, y12, 2, 1.e-13, 1.e-12);
	for(int i = 0; i < size; i++) printf("%lf ", (y8[i] - y10[i])/(y10[i] - y12[i]));
	printf("\n");
}



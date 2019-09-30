#include "h.h"
#define NUM_POINTS 100

void f(double x, double* y, int size, double* res) {
	res[1] = y[0] - y[1];
	res[0] = y[0] - y[1] - x/y[1];
}

void plot(char* fname, int num_points) {
	FILE* gnuplot_pipe = popen("gnuplot -persistent", "w");
	fprintf(gnuplot_pipe, "%s '%s' %s\n", "plot", fname, "with line");
}

int main(void) {
	double y0[2], y[2], a = 1., b = 0., h = (b - a)/NUM_POINTS, x = a;
	y0[0] = 0.; y0[1] = 0.;
        shooting(y0, 2, 1, a, b);
	y0[1] = 0.411;
	FILE* f = fopen("track.txt", "w");
	for(int i = 0; i < NUM_POINTS; i++, x += h) {
		runge(a, x, y0, y, 2, 1000);
		fprintf(f, "%lf %lf \n", y[0], y[1]);	
	}
	plot("track.txt", NUM_POINTS);
	fclose(f);
	return 0;	
}

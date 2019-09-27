#include "h.h"
#define NUM_POINTS 100

void plot(double* y, int num_points);
void f(double x, double* y, int size, double* res) {
	res[0] = y[1];
	res[1] = -1.*y[0];
}

void plot(double* y, int num_points) {
        char* gnuplot_commands[] = {"set title \"curve\"", "plot 'data.temp' with line"};
        FILE * temp = fopen("data.temp", "w");
        FILE * gnuplot_pipe = popen("gnuplot -persistent", "w");
        int i;
        for (i = 0; i < num_points; i++)
                fprintf(temp, "%lf %lf \n", y[2*i], y[2*i+1]);
        fprintf(gnuplot_pipe, "%s \n", gnuplot_commands[0]);
        fprintf(gnuplot_pipe, "%s \n", gnuplot_commands[1]);
}

int main(void) {
	double x0 = 0, x = 0, y0[2], y[2*NUM_POINTS], h = 2*3.14/NUM_POINTS;
	y0[0] = 0; y0[1] = 1;
	for(int i = 0; i < NUM_POINTS; i++, x+=h) {
		//euler(x0, x, y0, y, 2, 10000);
		runge_with_autostep(x0, x, y0, y + 2*i, 2,0.0000001, 0.00001);
	}
	plot(y, NUM_POINTS);
	return 0;
}



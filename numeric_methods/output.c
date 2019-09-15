#include "h.h"

void print(double* y, int size) {
	int i;
	for(i = 0; i < size; i++ )
		printf(" %lf ", y[i]);
	printf("\n");
}
int plot(double* y, int num_points) {
	char* gnuplot_commands[] = {"set title \"curve\"", "plot 'data.temp'"};
	FILE * temp = fopen("data.temp", "w");
	FILE * gnuplot_pipe = popen("gnuplot -persistent", "w");
	int i;
	for (i = 0; i < num_points; i++)
		fprintf(temp, "%lf %lf \n", y[2*i], y[2*i+1]);
	fprintf(gnuplot_pipe, "%s \n", gnuplot_commands[0]);
	fprintf(gnuplot_pipe, "%s \n", gnuplot_commands[1]);
	return 0;
}

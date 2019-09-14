#include "h.h"

void f(double x, double* y, int size, double* res) {
	int i;
	for(i = 0; i < size; i++) res[i] = x;
}

int main(void) {
	int i;
	double x0 = 0, x = 5, y0[1], y[1];
	y0[0] = 0; 
	evaluate_with_autostep(x0, x, y0, y, 1, 0.02, 0.000001, -1);	
	printf("y =  %lf\n", y[0]);
	return 0;
}

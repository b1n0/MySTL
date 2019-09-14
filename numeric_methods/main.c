#include "h.h"

void f(double x, double* y, int size, double* res) {
	int i;
	for(i = 0; i < size; i++) res[i] = x;
}

int main(void) {
	int i;
	double x0 = 0, x = 2, y0[1], y[1];
	y0[0] = 0; 
	for( i = 10; i < 25; i++) {
		evaluate(x0, x, y0, y, 1, i);	
		printf("steps = %d, y =  %lf\n", i, y[0]);
	}
	return 0;
}

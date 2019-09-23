#include "h.h"
#define NUM_POINTS 100

void f(double x, double* y, int size, double* res) {
	res[0] = y[1];
	res[1] = -1.*y[0];
}
		
int main(void) {
	int i;
	double x0 = 0, x = 0, y0[2], y[2*NUM_POINTS], h = 0.1;
	y0[0] = 1; y0[1] = 0; 
	for(i = 0; i < NUM_POINTS; i++, x+=h) 
		//euler(x0, x, y0, y + 2*i, 2, 10000);
		runge(x0, x, y0, y0 + 2*i, 2, 100);
		//runge_with_autostep(x0, x, y0, y + 2*i, 2, 0.02, 0.00000001, 10);
	plot(y, NUM_POINTS);	
	return 0;
}

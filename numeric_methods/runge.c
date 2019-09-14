#include "h.h"
#define sum(y, k, c) for(i = 0; i < size; i++) buff[i] = (y)[i] + (c)*(k)[i] 

int evaluate(double x0, double x, double* y0, double *y, int size, int num_steps) {
	int i = 0;
	double h = (x - x0)/num_steps, k[4*ST_SIZE], buff[ST_SIZE];
	memcpy(y, y0, size*sizeof(double)); 
	for (; x0 < x - h; x0 += h) {
		f(x0, y, size, k);
		sum(y, k, h/2);
		f(x0 + h/2, buff, size, k + size);
		sum(y, k + size, h/2);
		f(x0 + h/2, buff, size, k + 2*size);
		sum(y, k + 2*size, h);
		f(x0 + h, buff, size, k + 3*size);	
	
		for(i = 0; i < size; i++) 
			y[i] = y[i] + (k[i] + 2*k[size + i] + 2*k[2*size + i] + k[3*size + i])*h/6;
	}
	return 0;
}


#include "h.h"

#define sum(y, k, c) for(i = 0; i < size; i++) buff[i] = (y)[i] + (c)*(k)[i] 
#define ST_SIZE 32

int runge(double x0, double x, double* y0, double *y, int size, int num_steps) {
	int i = 0, j = 0;
	double h = (x - x0)/num_steps, k[4*ST_SIZE], buff[ST_SIZE];
	memcpy(y, y0, size*sizeof(double)); 
	for (j = 0; j < num_steps; j++) {
		f(x0, y, size, k);
		sum(y, k, h/2);
		f(x0 + h/2, buff, size, k + size);
		sum(y, k + size, h/2);
		f(x0 + h/2, buff, size, k + 2*size);
		sum(y, k + 2*size, h);
		x0 += h;
		f(x0, buff, size, k + 3*size);	
		for(i = 0; i < size; i++) 
			y[i] = y[i] + (k[i] + 2.*k[size + i] + 2.*k[2*size + i] + k[3*size + i])*h/6.;
	}
	return 0;
}

int runge_with_autostep(double x0, double x, double* y0, double* y, int size, double h, double err, double K) {
	int i = 0, j = 0;
	double k[4*ST_SIZE], buff[ST_SIZE], E;
	memcpy(y, y0, size * sizeof(double));
	for(; x0 < x - h; x0 += h) {
		while (1) {
			f(x0, y, size, k);
			sum(y, k, h/2);
			f(x0 + h/2, buff, size, k + size);
			sum(y, k + size, h/2);
			f(x0 + h/2, buff, size, k + 2*size);
			sum(y, k + 2*size, h);
			f(x0 + h, buff, size, k + 3*size);
			for(E = 0, i = 0; i < size; i++)
				E +=(k[i] - k[size + i] - k[2*size + i] + k[3*size + i]) * 
					(k[i] - k[size + i] - k[2*size + i] + k[3*size + i]) * 4./9.;
			if (E < err/K)  h *= 2;
			else if (E > err)  h /= 2; 
			else break;
		}
		for(i = 0; i < size; i++)
			y[i] += (k[i] + 2.*k[size + i] + 2.*k[2*size + i] + k[3*size + i])*h/6.;
	}
	runge(x0, x, y, buff, size, 20);
	memcpy(y, buff, size * sizeof(double));
	return 0;
}

int euler(double x0, double x, double* y0, double* y, int size, int num_steps) {
	int i, j;
	double h = (x - x0)/num_steps, buff[ST_SIZE];
	memcpy(y, y0, size * sizeof(double));
	for(i = 0; i < num_steps; i++, x0+=h) {
		f(x0, y, size, buff);
		for(j = 0; j < size; j++)
			y[j] += h*buff[j];
	}
	return 0;
}

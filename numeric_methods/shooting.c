#include "h.h"

#define delta 0.00000001
#define eps 0.0000001
#define random(min, max) (min) + ((double)rand()/RAND_MAX)*((max) - (min))

int shooting(double* y0, int size, int k, double a, double b) {
	double **m, *A, y0_buff[ST_SIZE], y[ST_SIZE], v[ST_SIZE], h[ST_SIZE], err;
	int i, j;

	m = (double**)malloc((size - k)*sizeof(double*));
	A = (double*)malloc((size - k)*size*sizeof(double));
	for(i = 0; i < size - k; i++) m[i] = m[i] = A + i*size + k;
	for(i = 0; i < size - k; i++) y0_buff[k+i] = random(-1., 1.);
	memcpy(y0_buff, y0, k*sizeof(double));
	
	while(1) {
		runge_with_autostep(a, b, y0_buff, y, size, 1.e-9, 1.e-8);
		for(i = 0 ; i < size - k; i++)
			v[i] = y0[k+i] - y[k+i];
		err = norm(v, size - k);
		if(err > eps) {
			for(i = 0; i < size - k; i++) {
				if (i > 0) 
		       			y0_buff[k + i - 1] -= delta; 
				y0_buff[k + i] += delta;
				runge_with_autostep(a, b, y0_buff, A + i*size, size, 1.e-9, 1.e-8);
				for(j = 0; j < size - k; j++) 
					m[i][j] = (m[i][j] - y[k+j])/delta;
			}
			y0_buff[size-1] -= delta;
			gauss(m, h, v, size - k);
			for(i = 0; i < size - k; i++) 
				y0_buff[k+i] += h[i];
		}
		else break;
	}
	free(A); free(m);
	memcpy(y0, y0_buff, size*sizeof(double));
	return 0;
}


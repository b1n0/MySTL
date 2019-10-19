#include "h.h"

int shooting(double a, double b, double* y0, int size, int k, double eps) {
	double **m, *A, y0_buff[ST_SIZE], y[ST_SIZE], v[ST_SIZE], h[ST_SIZE], err;
	int i, j;
	m = (double**)malloc((size - k)*sizeof(double*));
	A = (double*)malloc((size - k)*size*sizeof(double));
	for(i = 0; i < size - k; i++) m[i] = m[i] = A + i*size + k;
	for(i = 0; i < size - k; i++) y0_buff[k+i] = 1.;
	memcpy(y0_buff, y0, k*sizeof(double));
	
	while(1) {
		runge_hardcore(a, b, y0_buff, y, size, 1.e-7, 1.e-6);
		for(i = 0 ; i < size - k; i++)
			v[i] = y0[k+i] - y[k+i];
		err = norm(v, size - k, 'm');
		printf("%lf \n", err);
		if(err > eps) {
			for(i = 0; i < size - k; i++) {
				if (i > 0) 
		       			y0_buff[k + i - 1] -= DELTA; 
				y0_buff[k + i] += DELTA;
				runge_hardcore(a, b, y0_buff, A + i*size, size, 1.e-7, 1.e-6);
				for(j = 0; j < size - k; j++) 
					m[i][j] = (m[i][j] - y[k+j])/DELTA;
			}
			y0_buff[size-1] -= DELTA;
			gauss(m, h, v, size - k);
			for(i = 0; i < size - k; i++) y0_buff[k+i] += h[i];
		}
		else break;
	}
	free(A); free(m);
	memcpy(y0, y0_buff, size*sizeof(double));
	return 0;
}

int shoot(double a, double b, double* y0, int size, int k, double eps) {
	int i, j, n;
	double **m, y0_buff[ST_SIZE], v[ST_SIZE], y[ST_SIZE], h[ST_SIZE];
	double err, prev_err, c;
	prev_err = 10.;
	m = create_matrix(size - k, size - k);
	for(i = k; i < size; i++) y0[i] = 1.;
	//add loop for different beta
	//add return 0 condition
	runge_hardcore(a, b, y0, y, size, 1.e-8, 1.e-7);
	discrepancy(y0, y, v);
	for(err = norm(v, size - k, 'm'); err > eps; prev_err = err) {
		jacobian(m, y0, v);
		gauss(m, h, v, size - k);
		for(c = 1., i = 0; i < 50; i++, c*=0.5 ) {
			for(j = k; j < size; j++) y0_buff[j] = y0[j] - c*h[j - k];	
			runge_hardcore(a, b, y0_buff, y, size, 1.e-8, 1.e-7);
			discrepancy(y0, y, v);
			err = norm(v, size - k, 'm');
			printf("%lf \n", err);
			if(err < prev_err) break;
		}
		memcpy(y0 + k, y0_buff + k, sizeof(double)*(size - k));
	}	
	delete_matrix(m, size - k);
	return 1;
}

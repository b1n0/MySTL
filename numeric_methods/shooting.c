#include "h.h"

int shoot(double a, double b, double* y0, int size, int k, double eps) {
	int i, j, n, flag;
	double **m, y0_buff[ST_SIZE], v[ST_SIZE], y[ST_SIZE], h[ST_SIZE];
	double err, prev_err, c;
	m = create_matrix(size - k, size - k);
	for(i = k; i < size; i++) y0[i] = 3.;
	//add loop for different beta
	//add return 0 condition
	y0_buff[0] = y0[0];
	memcpy(y0_buff, y0, sizeof(double)*size);
	runge_hardcore(a, b, y0, y, size, 1.e-11, 1.e-10);
	discrepancy(y0, y, v);
	for(err = norm(v, size - k, 'm'), prev_err = err; err > eps; prev_err = err) {
		printf("%lf \n", err);
		for(i = k; i < size; i++) {
			y0[i] += DELTA;
			runge_hardcore(a, b, y0, y, size, 1.e-11, 1.e-10);
			discrepancy(y0, y, h);
			for(j = 0; j < size - k; j++) m[j][i-k] = (h[j] - v[j])/DELTA;
			y0[i] -= DELTA;
		}
		gauss(m, h, v, size - k);
		for(c = 1., flag = 1, i = 0; flag == 1 && i < 30 ; i++, c*=0.5 ) {
			for(j = k; j < size; j++) y0_buff[j] = y0[j] - c*h[j - k];	
			runge_hardcore(a, b, y0_buff, y, size, 1.e-11, 1.e-10);
			discrepancy(y0, y, v);
			err = norm(v, size - k, 'm');
			if(err < prev_err) flag = 0;
		}
		memcpy(y0 + k, y0_buff + k, sizeof(double)*(size - k));
	}	
	delete_matrix(m, size - k);
	return 1;
}

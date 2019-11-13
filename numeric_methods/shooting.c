#include "h.h"

int shoot(double a, double b, double* y0, int size, int k, double eps, 
			void discrepancy(double* y0, double* y, double* res), double err_min, double err_max) {
	int i, j, flag, num_c, res = 0;
	double **m, y0_buff[ST_SIZE], v[ST_SIZE], y[ST_SIZE], h[ST_SIZE];
	double err, prev_err, c;
	m = create_matrix(size - k, size - k);
	start_value(y0);
	memcpy(y0_buff, y0, size*sizeof(double));
	integrate_autostep(a, b, y0, y, size, err_min, err_max, 1.e-2, dormand8, 0);
	discrepancy(y0, y, v);
	for(c = 1., num_c = 0, err = norm(v, size - k, 'm'), prev_err = err; err > eps; prev_err = err) {
		printf("%lf \n", err);
		for(i = k; i < size; i++) {
			y0[i] += DELTA;
			integrate_autostep(a, b, y0, y, size, err_min, err_max, 1.e-2, dormand8, 0);
			discrepancy(y0, y, h);
		for(j = 0; j < size - k; j++) m[j][i-k] = (h[j] - v[j])/DELTA;
			y0[i] -= DELTA;
		}
		gauss(m, h, v, size - k);
		//for(c = MIN(1, 2*c), flag = 1, num_c -= num_c > 0 ? 1 : 0; flag == 1 && num_c < 50 ; c*=0.5, num_c++) {
		for(c = 1., flag = 1, num_c = 0; flag == 1 && num_c < 50 ; c*=0.5, num_c++) {
			for(j = k; j < size; j++) y0_buff[j] = y0[j] - c*h[j - k];	
			integrate_autostep(a, b, y0_buff, y, size, err_min, err_max, 1.e-2, dormand8, 0);
			discrepancy(y0, y, v);
			err = norm(v, size - k, 'm');
			if(err < prev_err) flag = 0;
		}
		if(flag == 1) { res = -1; break; }
		memcpy(y0 + k, y0_buff + k, sizeof(double)*(size - k));
	}	
	delete_matrix(m, size - k);
	return res;
}

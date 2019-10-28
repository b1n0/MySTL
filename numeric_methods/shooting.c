#include "h.h"

int shoot(double a, double b, double* y0, int size, int k, double eps, 
			void discrepancy(double* y0, double* y, double* res)) {
	int i, j, n, flag, num_c, res = 0;
	double **m, y0_buff[ST_SIZE], v[ST_SIZE], y[ST_SIZE], h[ST_SIZE];
	double err, prev_err, c;
	m = create_matrix(size - k, size - k);
	start_value(y0);
	memcpy(y0_buff, y0, size*sizeof(double));
	runge_hardcore(a, b, y0, y, size, 1.e-8, 1.e-7);
	discrepancy(y0, y, v);
	for(c = 1., num_c = 0, err = norm(v, size - k, 'm'), prev_err = err; err > eps; prev_err = err) {
		printf("%lf \n", err);
		for(i = k; i < size; i++) {
			y0[i] += DELTA;
			runge_hardcore(a, b, y0, y, size, 1.e-8, 1.e-7);
			discrepancy(y0, y, h);
		for(j = 0; j < size - k; j++) m[j][i-k] = (h[j] - v[j])/DELTA;
			y0[i] -= DELTA;
		}
		gauss(m, h, v, size - k);
		for(c = MIN(1, 2*c), flag = 1, num_c -= num_c > 0 ? 1 : 0; flag == 1 && num_c < 50 ; c*=0.5, num_c++) {
			for(j = k; j < size; j++) y0_buff[j] = y0[j] - c*h[j - k];	
			runge_hardcore(a, b, y0_buff, y, size, 1.e-8, 1.e-7);
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


double g(double a, double b, double *y0, int size, int k) {
	double y[ST_SIZE], v[ST_SIZE], res = 0.;
	runge_hardcore(a, b, y0, y, size, 1.e-8, 1.e-7);
	discrepancy(y0, y, v);
	for(int i = 0; i < size - k; i++) res += v[i]*v[i];
	return res;
}
int gradient_decrease(double a, double b, double* x, int size, int k, double eps) {
	int i = 0;
	double gradient[ST_SIZE], val, prev_val = 0.;
	start_value(x);
	for(val = g(a, b, x, size, k); fabs(val - prev_val) > eps || i == 0; prev_val = val, val = g(a, b, x, size, k)) {
		for(i = k; i < size; i++) {
			x[i] += DELTA;
			gradient[i-k] = (g(a, b, x, size, k) - val)/DELTA;
			x[i] -= DELTA;
		}
		for(i = k; i < size; i++) x[i] -= gradient[i-k];
		printf("%lf\n", val);
	}
	return 0;
}


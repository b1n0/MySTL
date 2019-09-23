#include "h.h"

#define eps 0.0000001
#define is_zero(a) ((a) < 0.0000001 && (a) > -0.0000001)

int gauss(double** m, double* x, double* b, int n); 
int gauss2(double* m, double* x, double* b, int n);

int gauss(double** m, double* x, double* b, int n) {
	int i, j, k, ind, c = 0, rank = 0;
	double* tmp;
	for(k = 0; k < n; k++) {
		for(ind = -1, i = c; i < n; i++)
			if (!is_zero(m[i][k])) { 
				ind = i;
				for(j = 0; j < k; j++)
					m[i][j] /= m[i][k];
				b[i] /= m[i][k];
			}
		for(i = c; i < n && ind != -1; i++) {
			if (!is_zero(m[i][k]) && i != ind) {
				for(j = k; j < n; j++)
					m[i][j] -= m[ind][j];
				b[i] -= b[ind];
			}
		}
		if (ind != -1 && ind != c) {
			tmp = m[ind]; m[ind] = m[c]; m[c] = tmp;
			b[ind] += b[c]; b[c] = b[ind] - b[c]; b[ind] = b[ind] - b[c];
		}
		c++;
	}

	for(k = 0; k < n; k++) {
		double s;
		if is_zero(m[n-1-k][n-1-k])  x[n-1-k] = 0;	
		else {
			for(s = 0., i = 0; i < k; i++)
				s-= x[n-1-i]*m[n-1-k][n-1-i];
			x[n-1-k] = (b[n-1-k] + s)/m[n-1-k][n-1-k];
		}
	}
	return 0;
}

int main(void) {
	double *m[2], b[2], x[2];
	m[0] = (double*)malloc(2*sizeof(double));
	m[1] = (double*)malloc(2*sizeof(double));
	b[0] = 1; b[1] = 1;
	m[0][0] = m[1][1] = 1; m[0][1] = m[1][0] = 1;;
	gauss2(m, x, b, 2);
	printf("%lf %lf\n", m[0][0], m[0][1]);
	printf("%lf %lf\n", m[1][0], m[1][1]);
	printf("%lf %lf\n", x[0], x[1]);
	printf("%lf %lf\n", b[0], b[1]);
	free(m[0]); free(m[1]);
	return 0;
}

int gauss2(double* m, double* x, double* b, int n) {
	int i, j, k, not_zero, rank = 0;
	for(k = 0; k < n; k++) {
		for(not_zero = -1, i = 0; i < n; i++)
			if (m[i*n+ k] > eps || m[i*n+ k] < -eps) {
				not_zero = i;	
				for(j = k; j < n; j++)
					m[i*n+ j] /= m[i*n+ k];
				b[i] /= m[i*n+ k];
			}
		x[k] = not_zero;
		for(i = 0; i < n; i++) {
			if (i != not_zero && (m[i*n+k] >  eps || m[i*n+k] < -eps)) {
				for(j = k; j < n; j++)
					m[i*n+ j] -= m[not_zero*n+ j];
				b[i] -= b[not_zero];
			}
		}
	}
	return 0;
}


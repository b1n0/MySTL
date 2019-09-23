#include "h.h"
#define is_zero(a) ((a) < 0.000000001 && (a) > -0.000000001)

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

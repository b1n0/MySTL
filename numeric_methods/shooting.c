#include "h.h"

#define delta 0.00000001
#define ST_SIZE 32
#define NUM_STEPS 1000

int fill_random(double* beta, int k) {
	for(int i = 0; i < k; i++)
		beta[i] = 0; //RANDOM
}

int shooting(double* y0, int size, int k, double a, double b) {
	double beta[ST_SIZE], y0_buff[ST_SIZE], res[ST_SIZE], **m;

	fill_random(beta, k);
	memcpy(y0_buff, y0, (size-k)*sizeof(double));
	memcpy(y0_buff + size - k, beta, k*sizeof(double));
	runge(a, b, y0_buff, res, size, NUM_STEPS);
	
	for(i = 0; i < k; i++) {
		if (i > 0)  y0_buff[n-k+i-1] -= delta; 
		y0_buff[n-k+i] += delta;
		runge(a, b, y0_buff, , size, NUM_STEPS);
		m[i] = A[i] + n - k;
		for(j = 0; j < k; j++) 
			A[i][n-k+i] = (A[i][n-k+j] - res[n-k+j])/delta;
	}
	buff[n-1] -= h;

	solve(m, x, , k);
	beta += sol;
	for(i = 0; i < size; i++) free(A[i]);
	return 0;
}


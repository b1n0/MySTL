#include "h.h"

#define delta 0.00000001
#define ST_SIZE 32
#define NUM_STEPS 1000

int fill_random(double* beta, int k) {
	for(int i = 0; i < k; i++)
		beta[i] = 0; //RANDOM
}

int shooting(double* y0, int size, int k, double a, double b) {
	double beta[ST_SIZE], buff[ST_SIZE], res[ST_SIZE], *A[ST_SIZE];
	for(i = 0; i < size; i++)
		A[i] = (double*)malloc(k*size*sizeof(double));
	fill_random(beta, k);
	memcpy(buff, y0, (size-k)*sizeof(double));
	memcpy(buff + size - k, beta, k*sizeof(double));
	runge(a, b, buff, res, size, NUM_STEPS);
	
	for(i = 0; i < k; i++) {
		if (i > 0)  buff[i-1] -= delta; 
		buff[i] += delta;
		runge(a, b, buff, A+i, size, NUM_STEPS);
		for(j = 0; j < k; j++) 
			A[size-k+j] = (A[size-k+j] - res[size-k+j])/delta;
	}
	beta[k-1] -= h;

	solve(A, x, b, k);
	beta += sol;
	for(i = 0; i < size; i++) free(A[i]);
	return 0;
}


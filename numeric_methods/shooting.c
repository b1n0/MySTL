#include "h.h"

#define delta 0.00000001
#define NUM_STEPS 1000
#define NUM_SHOTS 1
#define random(min, max) (min) + ((double)rand()/RAND_MAX)*((max) - (min))

void print_vector(double *y, int size);

void print_vector(double* y, int size) {\
	for(int i = 0; i < size; i++)
		printf("%lf ", y[i]);
	printf("\n");
}

int fill_random(double* beta, int k) {
	for(int i = 0; i < k; i++)
		beta[i] = random(-1, 1); 
}

int shooting(double* y0, int size, int k, double a, double b) {
	double **m, y0_buff[ST_SIZE], res[ST_SIZE], v[ST_SIZE], h[ST_SIZE];
	int i, j;

	m = (double**)malloc((size-k)*sizeof(double*));
	for(i = 0; i < size - k; i++) m[i] = malloc(size*sizeof(double));
	fill_random(y0_buff + k, size - k);
	memcpy(y0_buff, y0, k*sizeof(double));
	
	for(int shot = 0; shot < NUM_SHOTS; shot++) {
		//runge(a, b, y0_buff, res, size, NUM_STEPS);
		runge_with_autostep(a, b, y0_buff, res, size, 0.00000001, 0.0001);
		for(i = 0; i < size - k; i++) {
			if (i > 0) 
		       		y0_buff[k + i - 1] -= delta; 
			y0_buff[k + i] += delta;
			//runge(a, b, y0_buff, m[i], size, NUM_STEPS);
			runge_with_autostep(a, b, y0_buff, m[i], size, 0.00000001, 0.0001);
			for(j = 0; j < size - k; j++) 
				m[i][k+j] = (m[i][k+j] - res[k+j])/delta;
			v[i] = y0[k+i] - res[k+i];
		}
		//print_vector(res, size);
		print_vector(v, size - k);
		y0_buff[size-1] -= delta;
		gauss(m, h, v, size - k);
		for(i = 0; i < size - k; i++) 
			y0_buff[k+i] += h[i];
	}
	for(i = 0; i < size - k; i++) free(m[i]);
	free(m);
	return 0;
}


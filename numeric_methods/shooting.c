#include "h.h"

#define delta 0.00000001
#define NUM_STEPS 1000
#define NUM_SHOTS 1
#define NUM_POINTS 100
#define random(min, max) (min) + ((double)rand()/RAND_MAX)*((max) - (min))

int shooting(double* y0, int size, int k, double a, double b) {
	double **m, y0_buff[ST_SIZE], res[ST_SIZE], v[ST_SIZE], h[ST_SIZE], x;
	int i, j;

	m = create_matrix(size - k, size);
	for(i = 0; i < size - k; i++)  y0_buff[k+i] = 0.43;//random(-1., 1.);
	memcpy(y0_buff, y0, k*sizeof(double));
	
	//while(norm(v) > eps) {	
	for(int shot = 0; shot < NUM_SHOTS; shot++) {
		runge_with_autostep(a, b, y0_buff, res, size, 0.000001, 0.0001);
		for(i = 0; i < size - k; i++) {
			if (i > 0) 
		       		y0_buff[k + i - 1] -= delta; 
			y0_buff[k + i] += delta;
			runge_with_autostep(a, b, y0_buff, m[i], size, 0.000001, 0.0001);
			for(j = 0; j < size - k; j++) 
				m[i][k+j] = (m[i][k+j] - res[k+j])/delta;
			v[i] = y0[k+i] - res[k+i];
		}
		print_vector(res, size, stdout);
		print_vector(y0_buff, size, stdout);
		y0_buff[size-1] -= delta;
		gauss(m, h, v, size - k);
		for(i = 0; i < size - k; i++) 
			y0_buff[k+i] += h[i];
	}
	FILE *f = fopen("boundary.txt", "w");
	FILE *yf = fopen("y.txt", "w");
	FILE *xf = fopen("x.txt", "w");
	for(double step = (b-a)/NUM_POINTS, x = a, i = 0; i < NUM_POINTS; i++, x+=step) {
		runge(a, x, y0_buff, res, size, 10000);
		print_vector(res, size, f);
		fprintf(xf, "%lf %lf \n", x, res[0]);
		fprintf(yf, "%lf %lf \n", x, res[1]);
	}
	fclose(f);
	fclose(yf);
	fclose(xf);
	delete_matrix(m, size - k);
	return 0;
}


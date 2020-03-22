#include "h.h"

void derivative(int n, double* x, void f(double* x, double* y), double** res, double* buff) {
	int i,j;
	f(x, buff);
	for(i = 0; i < n-1; i++) {
		x[i] += DELTA;
		f(x, buff + n);
	       	x[i] -= DELTA;
		for(j = 0; j < n; j++) 
			res[j][i] = (buff[j + n] - buff[j])/DELTA;	
	}
}

int newton_method(int n, double* x, void f(double* x, double* y), double eps) {
	int i;
	double** jac = create_matrix(n, n), *buff = (double*)malloc(3*n*sizeof(double));
	f(x, buff + n);
	while(norm(buff + n, n, 'm') > eps) {
		jacobian(n, x, jac, buff + 2*n);
		gauss(jac, buff, buff + n, n);
		for(i = 0; i < n; i++) 
			x[i] += buff[i];
		f(x, buff + n);
	}		
	delete_matrix(jac, n);
	free(buff);
	return 1;
}

void f(double* x, double* y) {
	y[0] = x[0] - 5.;
}

int main(void) {
	double* x = (double*)malloc(sizeof(double));
	newton_method(1, x, f, 1.e-5);
	return 0;
}

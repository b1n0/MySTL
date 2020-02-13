#include "h.h"
#define D(x, u, du, d2u) (a*(u) + 1. + (du)/(x) + (d2u))

void grid_evaluate(double** u, int tn, int xn, double a, double b) {
	double x, t, h = 1./(xn-1), tau = 1./(tn-1);
	int i, j;
	for(i = xn - 1., x = 1.; i > 0; i--, x -= xstep) 
		u[0][i] = tau*D(x, 0.5*(1-x*x), -x, -1.) + 0.5*(1-x*x) 
	u[0][0] = h*h*(u[0][1] - 0.5*(1-h*h))/tau - u[0][2] + 2.*u[0][1];
	for(i = 1; i < tn; i++) {
		for(j = xn - 1, x = 1.; j > 0; j--, x -= xstep) 
			u[i][j] = 
		u[i][0] = h*h*(u[i][1] - u[i-1][1])/tau - u[i][2] + 2.*u[i][1]
	}
}

int main(void) {
	int xn, tn, i;
	double a, b, **u;
	FILE* f = fopen("out.txt", "w");

	u = create_matrix(tn, xn);
	grid_evaluate(u, tn, xn, a, b);
	print(u);

	delete_matrix(u);
	fclose(f);
	return 0;
}

#include "h.h"

void explicit_solve(double** u, int tn, int xn, double a) {
	double x, t, h = 1./(xn-1), tau = 1./(tn-1);
	int i, j;
	for(i = 0, x = 0.; i < xn; i++, x += h) u[0][i] = 0.5*(1. - x*x);		
	for(i = 0; i < tn-1; i++) {
		u[i+1][xn-1] = 0;
		u[i+1][0] = tau*(a*u[i][0] + 1.) + u[i][0];
		for(j = 1, x = h; j < xn - 1; j++, x += h)
			u[i+1][j] = tau*((u[i][j+1] - u[i][j-1])/(2.*x*h) + (u[i][j+1] -2.*u[i][j] + u[i][j-1])/(h*h) + a*u[i][j] + 1.) + u[i][j];
	}
}

void implicit_solve(double** u, int tn, int xn, double a) {
	double x, h = 1./(xn-1), tau = 1./(tn-1);
	double** A = create_matrix(xn - 2, xn - 2), *b = (double*)malloc((xn-2)*sizeof(double));
	int i,j;
	for(i = 0, x = 0.; i < xn; i++, x += h) u[0][i] = 0.5*(1. - x*x);
	for(i = 0; i < tn-1; i++) {
		u[i+1][xn-1] = 0.;
		u[i+1][0] = tau*(a*u[i][0] + 1.) + u[i][0];
		for(j = 0, x = h; j < xn - 2; j++, x+=h) {
			if(j != 0) A[j][j-1] = 1./(h*h) - 1./(2.*h*x); 
			if(j != xn - 3) A[j][j+1] = 1./(h*h*x) + 1./(2.*x*h);
			A[j][j] = (a - 2./(h*h) - 1./tau);
			b[j] = -u[i][j]/tau -1.;
		}
		gauss(A, u[i+1] + 1, b, xn - 2);
	}
	delete_matrix(A, xn-2); free(b);
}

int main(void) {
	int xn, tn, type;
	double a, **u;
	FILE* f = fopen("out.txt", "w");

	printf("enter 1 for explicit and 2 for implicit. xn tn alpha.\n");
	scanf("%d %d %d %lf", &type, &xn, &tn, &a);
	u = create_matrix(tn, xn);

	if(type == 1) explicit_solve(u, tn, xn, a);
	else implicit_solve(u, tn, xn, a);

	for(int i = tn - 1; i >= 0; i--) {
		for(int j = 0; j < xn; j++) fprintf(f, "%lf ", u[i][j]);
		fprintf(f, "\n");
	}

	delete_matrix(u, tn);
	fclose(f);
	return 0;
}

#include "h.h"

#define Tn 10001
#define Xn 10001

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

void running_method(int n, double* x, double** c, double* f) {
	double *a, *b;
	a = (double*)malloc(n*sizeof(double));
	b = (double*)malloc(n*sizeof(double));
	a[0] = c[1][0]/c[2][0];
	b[0] = f[0]/c[2][0];
	for(int i = 0; i < n - 1; i++) {
		a[i+1] = c[1][i]/(c[2][i] - c[0][i]*a[i]);
		b[i+1] = (f[i] + c[0][i]*b[i])/(c[2][i] - c[0][i]*a[i]);
	}
	x[n-1] = (f[n-1] + c[0][n-1]*b[n-1])/(c[2][n-1] - c[0][n-1]*a[n-1]);
	for(int i = n - 2; i >= 0; i--) 
		x[i] = a[i+1]*x[i+1] + b[i+1];
	free(a); free(b);
}

void implicit_solve(double** u, int tn, int xn, double a) {
	int i,j, n = xn - 1;
	double x, h = 1./(xn-1), tau = 1./(tn-1);
	double **c = create_matrix(3, n), *f = (double*)malloc(n*sizeof(double));

	c[0][0] = 0.; c[1][0] = 1.; c[2][0] = 1.; f[0] = 0.;
	for(j = 1, x = h; j < xn - 1; j++, x+=h) {
		c[0][j] = -xn + 1.5 + 0.5/x; 
		c[1][j] = -xn + 0.5 - 0.5/x;
		c[2][j] = a*h - 2.*xn + 2. - h*(tn - 1.);
	}
	for(j = 0, x = 0.; j < xn; j++, x += h) u[0][j] = 0.5*(1. - x*x);
	for(i = 0; i < tn-1; i++) {
		u[i+1][xn-1] = 0.;
		for(j = 1, x = h; j < xn - 1; j++, x+=h) f[j] = -h*tn*u[i][j];
	running_method(n, u[i+1], c, f);
	}
	free(f); delete_matrix(c, 3);
}

void inaccuracy(double **u, double** v, int tn, int xn, double a, double* ans) {
	int i, j, k, l, tstep = (int)((Tn-1)/(tn-1)), xstep = (int)((Xn-1)/(xn-1));
	double mnorm = 0., l1norm = 0., vl1norm = 0., vmnorm = 0.; 
	for(k = i = 0; i < tn; i++, k += tstep)
		for(l = j = 0; j < xn; j++, l += xstep) {
			l1norm += fabs(u[i][j] - v[k][l]);
			mnorm = MAX(mnorm, fabs(u[i][j] - v[k][l])); 
			vl1norm += fabs(v[k][l]);
			vmnorm = MAX(vmnorm, v[k][l]);
		}
	ans[0] = mnorm; ans[1] = l1norm/(xn-1); 
	ans[2] = mnorm/vmnorm; ans[3] = l1norm/vl1norm;
}

void save(double** u, int tn, int xn, const char* fname) {
	FILE* f = fopen(fname, "w");
	for(int i = tn - 1; i >= 0; i--) {
		for(int j = 0; j < xn; j++) fprintf(f, "%lf ", u[i][j]);
		fprintf(f, "\n"); 
	}
	fclose(f);
}

int main(void) {
	int i, n;
	double ans[4], **u, **v;
	v = create_matrix(Tn, Xn);
	implicit_solve(v, Tn, Xn, 1.);

	for(i = 0, n = 11; i < 3; i++, n = (n-1)*10 + 1) {
		u = create_matrix(n, n);
		implicit_solve(u, n, n, 1.);
		inaccuracy(u, v, n, n, 1., ans);
		printf("%lf %lf %lf %lf \n", ans[0], ans[1], ans[2], ans[3]);
		delete_matrix(u, n);
	}
	delete_matrix(v, Tn);
	return 0;
}

/*
int main(void) {
	int xn, tn, type;
	double a, ans[4], **u;
	double **v = create_matrix(Tn, Xn);

	printf("enter 1 for explicit and 2 for implicit, tn, xn, alpha.\n");
	scanf("%d %d %d %lf", &type, &tn, &xn, &a);
	u = create_matrix(tn, xn);
	
	if(type == 1) explicit_solve(u, tn, xn, a);
	else implicit_solve(u, tn, xn, a);
	implicit_solve(v, Tn, Xn, a);
	inaccuracy(u, v, tn, xn, a, ans);
	printf("%lf %lf %lf %lf \n", ans[0], ans[1], ans[2], ans[3]);
	save(u, tn, xn, "out.txt");

	delete_matrix(u, tn);
	return 0;
}
*/


#include "grid.h"

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
	int i,j, n = xn - 1;
	double x, h = 1./(xn-1), tau = 1./(tn-1);
	double **c = create_matrix(3, n), *f = (double*)malloc(n*sizeof(double));

	c[0][0] = 0.; c[1][0] = -1.; c[2][0] = -1.; f[0] = 0.;
	for(j = 1, x = h; j < xn - 1; j++, x+=h) {
		c[0][j] = -xn + 1.5 + 0.5/x;
		c[1][j] = j != xn - 2 ? -xn + 0.5 - 0.5/x : 0.;
		c[2][j] =  a*h - 2.*xn + 2. - h*(tn - 1.);
	}
	for(j = 0, x = 0.; j < xn; j++, x += h) u[0][j] = 0.5*(1. - x*x);
	for(i = 0; i < tn-1; i++) {
		u[i+1][xn-1] = 0.;
		for(j = 1, x = h; j < xn - 1; j++, x+=h) f[j] = -h*tn*u[i][j];
		running_method(n, u[i+1], c, f);
	}
	free(f); delete_matrix(c, 3);
}

void inaccuracy(double* u, double* v, int xn, double a, double* ans) {
	int i, k, xstep = (int)((Xn-1)/(xn-1));
	double max, l1, vl1, vmax;
	max = l1 = vl1 = vmax = 0.; 
	for(k = i = 0; i < xn; i++, k += xstep) {
		l1 += fabs(u[i] - v[k]);
		max = MAX(max, fabs(u[i] - v[k]));
		vl1 += fabs(v[k]);
		vmax = MAX(vmax, v[k]);
	}
	ans[0] = max; ans[1] = l1/(xn-1); 
	ans[2] = max/vmax; ans[3] = l1/vl1;
}

void save(double** u, int tn, int xn, const char* fname) {
	FILE* f = fopen(fname, "w");
	for(int i = tn - 1; i >= 0; i--) {
		for(int j = 0; j < xn; j++) fprintf(f, "%lf ", u[i][j]);
		fprintf(f, "\n"); 
	}
	fclose(f);
}

void plot(double** u, int tn, int xn, const char* fname) {
	int i,j;
	double tau = 1./(tn-1), h = 1./(xn-1), t, x, f1;
	FILE *gnuplot_pipe, *f = fopen(fname, "w");
	for(i = 0, t = 0.; i < tn; i++, t += tau) {
		for(j = 0, x = 0, f1 = 0.; j < xn; j++, f1 += u[i][j]*x*h, x += h);
		fprintf(f, "%lf %lf %lf\n", t, f1, (u[i][xn-1] - u[i][xn-2])*(xn-1));
	}
	fclose(f);
	gnuplot_pipe = popen("gnuplot -persistent", "w");
        fprintf(gnuplot_pipe, "%s%s%s\n", "set key off; plot for[col=2:3:1]'", fname, "' using 1:col with lines");
	fclose(gnuplot_pipe);

	f = fopen(fname, "w");
	gnuplot_pipe = popen("gnuplot -persistent", "w");
	for(i = 0, x = 0.; i < Xn; i++, x += h) 
		fprintf(f, "%lf %lf %lf %lf %lf\n", x, u[(int)Tn/4][i], u[(int)Tn/2][i], u[(int)3*Tn/4][i], u[Tn-1][i]);
	fclose(f);
        fprintf(gnuplot_pipe, "%s%s%s\n", "set key off; plot for[col=2:5:1]'", fname, "' using 1:col with lines");
	fclose(gnuplot_pipe);
}


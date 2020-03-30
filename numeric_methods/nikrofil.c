#include "h.h"
#define D(u, d2u) (alpha*(d2u) + pow((u), 4))

double alpha = 1.;
double beta = 1.;

void explicit_evaluate(double** u, int tn, int xn) {
	int i,j;
	double x, xstep = 1./(xn-1), tstep = 1./(tn-1);
	for(i = 0, x = 0.; i < xn; i++, x += xstep) 
		u[0][i] = beta*(1. - x*x)*(1. - x*x);
	for(i = 0; i < tn - 1; i++) {
		u[i+1][0] = tstep*D(u[i][0], (u[i][1] - u[i][0])/(xstep*xstep)) + u[i][0];
		u[i+1][xn-1] = tstep*D(u[i][xn-1], (-50.*pow(u[i][xn-1], 4) - (u[i][xn-1] - u[i][xn-2])/xstep)/xstep) + u[i][xn-1];   
		for(j = 1, x = 0.; j < xn-1; j++, x += xstep)
			u[i+1][j] = tstep*D(u[i][j], (u[i][j+1] - 2*u[i][j] + u[i][j-1])/(xstep*xstep)) + u[i][j]; 
	}
}

void f(int n, double* x, double* y, double* params, double h, double tau) {
	for(int i = 0; i < n; i++) 
		y[i] = alpha*x[i+1]/(h*h) - (2.*alpha/(h*h) + 1./tau)*x[i] + alpha*x[i-1]/(h*h) + params[i] + pow(x[i], 4);		
}

void jacobian(int n, double* x, double** res, double h, double tau) {
	int i, j;
	for(i = 0; i < n; i++) {
		if(i != 0) res [i][i-1] = alpha/(h*h);
		if(i != n-1) res[i][i+1] = alpha/(h*h);
		res[i][i] = 2.*alpha/(h*h) + 1./tau + 4.*pow(x[i], 3.);
	}
}

int implicit_evaluate(double** u, int tn, int xn, double eps) {
	int i,j,k, flag, n = xn - 2, num_c;
	double x, xstep = 1./(xn-1), tstep = 1./(tn-1), err, prev_err, c;
	double** jac = create_matrix(n, n), **v = create_matrix(3, n);
	for(i = 0, x = 0.; i < xn; i++, x += xstep) 
		u[0][i] = beta*(1. - x*x)*(1. - x*x);
	for(i = 0; i < tn - 1; i++) {
		u[i+1][0] = tstep*D(u[i][0], (u[i][1] - u[i][0])/(xstep*xstep)) + u[i][0];
		u[i+1][xn-1] = tstep*D(u[i][xn-1], (-50.*pow(u[i][xn-1], 4) - (u[i][xn-1] - u[i][xn-2])/xstep)/xstep) + u[i][xn-1];
		f(n, u[i+1] + 1, v[0], u[i] + 1, xstep, tstep);
		for(err = norm(v[0], n, 'm'), prev_err = err; err > eps; prev_err = err) {
			jacobian(n, u[i+1] + 1, jac, xstep, tstep);
			gauss(jac, v[1], v[0], n);
			for(flag = 1, c = 1., num_c = 0; flag && num_c < 50; num_c++, c *= 0.5) {
				for(k = 0; k < n; k++) v[2][k] = u[i+1][k+1] - c*v[1][k];
				f(n, v[2], v[0], u[i] + 1, xstep, tstep);
				err = norm(v[0], n, 'm');
				if(err < prev_err) flag = 0;
			}
			if(flag) return -1;
			memcpy(u[i+1] + 1, v[2], n*sizeof(double));
		} 
	}
	delete_matrix(jac, n); delete_matrix(v, 3);
}

int main(void) {
	double **u;
	int i,j, type, tn = 10, xn = 10;
	FILE * f = fopen("out.txt", "w");
	u = create_matrix(tn, xn);
	
	scanf("%d", &type);
	if(type == 1) explicit_evaluate(u, tn, xn);
	else implicit_evaluate(u, tn, xn, 1.e-5);
	for(i = tn - 1; i >= 0; i--) {
		for(j = 0; j < xn; j++) fprintf(f, "%lf ", u[i][j]);
		fprintf(f, "\n");
	}
	delete_matrix(u, tn); fclose(f);
	return 0;
}


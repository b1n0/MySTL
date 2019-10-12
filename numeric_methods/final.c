#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"

#define EPS 0.000001
#define DELTA 1.e-8
#define MIN_STEP 0.000001
#define ST_SIZE 32
#define NUM_POINTS 300
#define MAX(a, b) (a)>(b)?(a):(b)

void plot(const char* fname);
void f(double x, double* y, int size, double* res);
void runge_numbers(double x0, double x, double* y0, int size);
double eigen_value(double x, double* y);
double runge(double x0, double x, double* y0, double *y, int size, int num_steps);
double runge_hardcore(double x0, double x, double* y0, double* y, int size, double err_min, double err_max);
int shoot(double a, double b, double* y0, double* y, int n);

void plot(const char* fname) {
	FILE* gnu_pipe = popen("gnuplot -persistent", "w");
	fprintf(gnu_pipe,"%s '%s' %s\n", "plot", fname, "with line");
}

void f(double x, double* y, int size, double* res) {
	size++;
	res[0] = y[1];
	res[1] = x/y[0];
}

double eigen_value(double x, double* y) { return 1 - x/(y[0]*y[0]); }

int shoot(double a, double b, double* y0, double* y, int n) {
	double u[2], err, d, prev_err, c = 1;
	for(int i = 0; i < n; i++) {
		runge_hardcore(a, b, y0, y, 2, 1.e-6, 1.e-5);
		//runge(a, b, y0, y, 2., 20000);
		err = y[0];
		if(fabs(err) < EPS) return 1;
		y0[1] += DELTA;
		runge_hardcore(a, b, y0, u, 2, 1.e-6, 1.e-5);
		//runge(a, b, y0, u, 2., 20000);
		d = (u[0] - err)/DELTA;
		y0[1] -= DELTA;
		if(i && fabs(err) > fabs(prev_err)) 
			c *= 0.5;
		y0[1] -= c*err/d;
		prev_err = err;
	}
	return 0;
}

int main(void) {
	int i, res;
        double y0[2], y[2], betas[10], a, b, lambda, beta, lambda_step, c, h, err = 0.;
	FILE *f = fopen("track.txt", "w");
	for(i = 0; i < 10; i++) betas[i] = 2*i + 1; 
        a = 0.; b = 1.;
        lambda = 0.0001;
        lambda_step = 0.1;
        c = 1.;
        while(lambda_step > MIN_STEP) {
		printf("%lf \n", lambda);
		for(y0[0] = lambda, y0[1] = betas[0], res = 0, i = 0; i < 10; i++, y0[0] = lambda, y0[1] = betas[i])
			res += shoot(b, a, y0, y, 100);
		if(res == 0) { lambda_step *= c; lambda += lambda_step; }
                else { c = 0.5; lambda_step *= c; lambda -= c*lambda_step; }
	}
	printf("min lambda = %5.30lf\n", lambda);
	y0[0] = lambda;
	shoot(b, a, y0, y, 100);
	printf("beta = %5.30lf\n", y0[1]);
	h = (b-a)/NUM_POINTS;
	for(double x = b; x > a; x -= h) {
		err += runge_hardcore(x, x-h, y0, y, 2, 1.e-7, 1.e-8);
		y0[0] = y[0]; y0[1] = y[1];
		fprintf(f, "%lf %lf\n", y[0], y[1]);
	}
	y0[0] = lambda; y0[1] = beta;
	printf("global error = %lf\n", err);
	plot("track.txt");
	
	runge_numbers(b, 0.75, y0, 2);
	runge_numbers(b, 0.5, y0, 2);
	runge_numbers(b, 0.25, y0, 2);
	runge_numbers(b, 0., y0, 2);
	
	fclose(f);
        return 0;
}

void runge_numbers(double x0, double x, double* y0, int size) {
        double y8[ST_SIZE], y10[ST_SIZE], y12[ST_SIZE];
        runge_hardcore(x0, x, y0, y8, size, 1.e-9, 1.e-8);
        runge_hardcore(x0, x, y0, y10, size, 1.e-11, 1.e-10);
        runge_hardcore(x0, x, y0, y12, size, 1.e-13, 1.e-12);

        printf("runge number in %lf \t ", x);
        for(int i = 0; i < size; i++) printf("%lf ", (y8[i] - y10[i])/(y10[i] - y12[i]));
        printf("\n");
}

double runge(double x0, double x, double* y0, double *y, int size, int num_steps) {
        int i = 0, j = 0;
        double h = (x - x0)/num_steps, k[6][ST_SIZE], buff[ST_SIZE], c, E, global_err= 0.;
        memcpy(y, y0, size*sizeof(double));
        for (j = 0; j < num_steps; j++, x0+=h) {
                f(x0, y, size, k[0]);
                for(i = 0; i < size; i++) buff[i] = y[i] + h*k[0][i]*0.5;
                f(x0 + h*0.5, buff, size, k[1]);
                for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i] + k[1][i])*0.25;
                f(x0 + h*0.5, buff, size, k[2]);
                for(i = 0; i < size; i++) buff[i] = y[i] + h*(2*k[2][i] - k[1][i]);
                f(x0 + h, buff, size, k[3]);
                for(i = 0; i < size; i++) buff[i] = y[i] + h*(7*k[0][i] +10*k[1][i] + k[3][i])/27;
                f(x0 + 2*h/3, buff, size, k[4]);
                for(i = 0; i < size; i++) buff[i] = y[i] + h*(28*k[0][i] - 125*k[1][i] + 546*k[2][i] + 54*k[3][i] - 378*k[4][i])/625;
                f(x0 + h*0.2, buff, size, k[5]);
                for(E = 0, i = 0; i < size; i++) {
                        c = ((-42)*k[0][i] - 244*k[2][i] - 21*k[3][i] + 162*k[4][i] + 125*k[5][i])*h/336;
			//E += fabs(c);
			E = MAX(fabs(E), fabs(c));
                }
                global_err = exp(h*eigen_value(x, y))*global_err + E;
                for(i = 0; i < size; i++)
                        y[i] += (14*k[0][i] + 35*k[3][i] + 162*k[4][i] + 125*k[5][i])*h/336;
        }
        return global_err;
}


double runge_hardcore(double x0, double x, double* y0, double* y, int size, double err_min, double err_max) {
	double k[7][ST_SIZE], y1[ST_SIZE], y2[ST_SIZE], buff[ST_SIZE], E, c, h = (x - x0)/1000000, global_err = 0.;
	int i;	
	memcpy(y, y0, size * sizeof(double));
	for(; (x0 < x - h && h > 0) || (x0 > x - h && h < 0); x0 += h) {
		while (1) {
			f(x0, y, size, k[0]);
			for(i = 0; i < size; i++) buff[i] = y[i] + h*k[0][i]*0.5;
			f(x0 + h*0.5, buff, size, k[1]);
			for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i] + k[1][i])*0.25;
			f(x0 + h*0.5, buff, size, k[2]);
			for(i = 0; i < size; i++) buff[i] = y[i] + h*(2*k[2][i] - k[1][i]);
			f(x0 + h, buff, size, k[3]);	
			for(i = 0; i < size; i++) buff[i] = y[i] + h*(7*k[0][i] +10*k[1][i] + k[3][i])/27;
			f(x0 + 2.*h/3., buff, size, k[4]);
			for(i = 0; i < size; i++) buff[i] = y[i] + h*(28*k[0][i] - 125*k[1][i] + 546*k[2][i] + 54*k[3][i] - 378*k[4][i])/625;
			f(x0 + h*0.2, buff, size, k[5]);
			
			for(i = 0; i < size; i++) y1[i] = y[i] + (14*k[0][i] + 35*k[3][i] + 162*k[4][i] + 125*k[5][i])*h/336;

			for(i = 0; i < size; i++) buff[i] = y[i] + h*k[0][i]*2./9. + h*k[1][i]*4./9.;
			f(x0 + h*2./3., y, size, k[2]);
			for(i = 0; i < size; i++) buff[i] = y[i] + h*k[0][i]*7./36. + h*k[1][i]*2./9. - h*k[2][i]/12.;
			f(x0 + h*1./3., y, size, k[3]);
			for(i = 0; i < size; i++) buff[i] = y[i] - h*k[0][i]*35./144. - h*k[1][i]*55./36. + h*k[2][i]*35./48. + h*k[3][i]*15./8.;
			f(x0 + h*5./6., y, size, k[4]);
			for(i = 0; i < size; i++) buff[i] = y[i] - h*k[0][i]/360. - h*k[1][i]*11./36. - h*k[2][i]*0.125 + h*k[3][i]*0.5 + h*k[4][i]*0.1;
			f(x0 + h/6., y, size, k[5]);
			for(i=0;i<size;i++) buff[i]=y[i]-h*k[0][i]*41./260.+h*k[1][i]*22./13.+h*k[2][i]*43./156.-h*k[3][i]*118./39.+h*k[4][i]*32./195.+h*k[5][i]*80./39.;
			f(x0 + h, y, size, k[6]);

			for(i = 0; i < size; i++) y2[i] = y[i]+(k[0][i]*13./200. + k[2][i]*11./40. + k[3][i]*11./40. + k[4][i]*0.16 + k[5][i]*0.16 + k[6][i]*13./200.)*h;
			
			for(i = 0, E = 0.; i < size; i++)  E = MAX(fabs(E), fabs(y1[i] - y2[i]));

			if (E < err_min) { global_err = exp(h*eigen_value(x, y))*global_err + E; h*=2.; break; }
			else if(E > err_max) h *= 0.5; 
			else { global_err = exp(h*eigen_value(x, y))*global_err + E; break; }
		}
		memcpy(y, y1, size*sizeof(double));
	}
	global_err += runge(x0, x, y, y2, size, 200);
	memcpy(y, y2, size*sizeof(double));
	return global_err;
}

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"

#define EPS 0.0001
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
double runge_with_autostep(double x0, double x, double* y0, double* y, int size, double err_min, double err_max);
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
	double u[2], err, d;
	for(int i = 0; i < n; i++) {
		runge(a, b, y0, y, 2, 10000);
		err = y[0];
		if(err*err < EPS) return 1;
		y0[1] += DELTA;
		runge(a, b, y0, u, 2, 10000);
		d = (u[0] - err)/DELTA;
		y0[1] -= DELTA;
		y0[1] += -err/d;
	}
	return 0;
}

int main(void) {
        double y0[2], y[2], a, b, lambda, beta, lambda_step, c, err, h;
	FILE *f = fopen("track.txt", "w");
	//system("wget https://raw.githubusercontent.com/b1n0/study/master/wolfram/tmp.png && clear");
	//system("gsettings set org.gnome.desktop.background picture-uri file://$PWD/tmp.png");
        a = 0.; b = 1.;
        lambda = 0.1;
        beta = 2.;
        lambda_step = 1;
        c = 1.;
        while(lambda_step > MIN_STEP) {
                y0[0] = lambda;
        	y0[1] = beta;
		printf("%lf \n", lambda);
		if(shoot(b, a, y0, y, 3000) == 0) { lambda_step *= c; lambda += lambda_step; }
                else { c = 0.5; lambda_step *= c; lambda -= c*lambda_step; }
	}
	printf("min lambda = %lf\n", lambda);
	printf("beta = %lf\n", y0[1]);
	h = (b-a)/NUM_POINTS;
	for(double x = b; x > a + h; x -= h) {
		err = runge_with_autostep(x, x-h, y0, y, 2, 1.e-8, 1.e-7);
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
        runge_with_autostep(x0, x, y0, y8, size, 1.e-9, 1.e-8);
        runge_with_autostep(x0, x, y0, y10, size, 1.e-11, 1.e-10);
        runge_with_autostep(x0, x, y0, y12, size, 1.e-13, 1.e-12);

        printf("runge number in %lf \t ", x);
        for(int i = 0; i < size; i++) printf("%lf ", (y8[i] - y10[i])/(y10[i] - y12[i]));
        printf("\n");
}

double runge_with_autostep(double x0, double x, double* y0, double* y, int size, double err_min, double err_max) {
        int i = 0;
        double k[6*ST_SIZE], buff[ST_SIZE], E, c, h = (x - x0)/100000, global_err = 0.;
        memcpy(y, y0, size * sizeof(double));
        for(; (x0 < x - h && h > 0) || (x0 > x - h && h < 0); x0 += h) {
                while (1) {
                        f(x0, y, size, k);
                        for(i = 0; i < size; i++) buff[i] = y[i] + h*k[i]*0.5;
                        f(x0 + h*0.5, buff, size, k + size);
                        for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[i] + k[size+i])*0.25;
                        f(x0 + h*0.5, buff, size, k + 2*size);
                        for(i = 0; i < size; i++) buff[i] = y[i] + h*(2*k[2*size+i] - k[size+i]);
                        f(x0 + h, buff, size, k + 3*size);
                        for(i = 0; i < size; i++) buff[i] = y[i] + h*(7*k[i] +10*k[size + i] + k[3*size + i])/27;
                        f(x0 + 2*h/3, buff, size, k + 4*size);
                        for(i = 0; i < size; i++) buff[i] = y[i] + h*(28*k[i] - 125*k[size+i] + 546*k[2*size+i] + 54*k[3*size+i] - 378*k[4*size+i])/625;
                        f(x0 + h*0.2, buff, size, k + 5*size);

                        for(E = 0, i = 0; i < size; i++) {
                                c = ((-42)*k[i] - 244*k[2*size+i] - 21*k[3*size+i] + 162*k[4*size+i] + 125*k[5*size+i])*h/336;
				E += fabs(c);
				//E = MAX(fabs(E), fabs(c));
                        }
                        if (E < err_min) { global_err = exp(h*eigen_value(x, y))*global_err + E; h*=2.; break;}
                        else if(E > err_max) h *= 0.5;
                        else { global_err = exp(h*eigen_value(x, y))*global_err + E; break; }
                }
                for(i = 0; i < size; i++)
                        y[i] += (14*k[i] + 35*k[3*size+i] + 162*k[4*size+i] + 125*k[5*size+i])*h/336;
        }
        global_err += runge(x0, x, y, buff, size, 200);
        memcpy(y, buff, size * sizeof(double));
        return global_err;
}

double runge(double x0, double x, double* y0, double *y, int size, int num_steps) {
        int i = 0, j = 0;
        double h = (x - x0)/num_steps, k[6*ST_SIZE], buff[ST_SIZE], c, E, global_err= 0.;
        memcpy(y, y0, size*sizeof(double));
        for (j = 0; j < num_steps; j++, x0+=h) {
                f(x0, y, size, k);
                for(i = 0; i < size; i++) buff[i] = y[i] + h*k[i]*0.5;
                f(x0 + h*0.5, buff, size, k + size);
                for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[i] + k[size+i])*0.25;
                f(x0 + h*0.5, buff, size, k + 2*size);
                for(i = 0; i < size; i++) buff[i] = y[i] + h*(2*k[2*size+i] - k[size+i]);
                f(x0 + h, buff, size, k + 3*size);
                for(i = 0; i < size; i++) buff[i] = y[i] + h*(7*k[i] +10*k[size + i] + k[3*size + i])/27;
                f(x0 + 2*h/3, buff, size, k + 4*size);
                for(i = 0; i < size; i++) buff[i] = y[i] + h*(28*k[i] - 125*k[size+i] + 546*k[2*size+i] + 54*k[3*size+i] - 378*k[4*size+i])/625;
                f(x0 + h*0.2, buff, size, k + 5*size);
                for(E = 0, i = 0; i < size; i++) {
                        c = ((-42)*k[i] - 244*k[2*size+i] - 21*k[3*size+i] + 162*k[4*size+i] + 125*k[5*size+i])*h/336;
			E += fabs(c);
			//E = MAX(fabs(E), fabs(c));
                }
                global_err = exp(h*eigen_value(x, y))*global_err + E;
                for(i = 0; i < size; i++)
                        y[i] += (14*k[i] + 35*k[3*size+i] + 162*k[4*size+i] + 125*k[5*size+i])*h/336;
        }
        return global_err;
}


#include "h.h"

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

double runge_with_autostep(double x0, double x, double* y0, double* y, int size, double err_min, double err_max) {
	int i = 0;
	double k[6*ST_SIZE], buff[ST_SIZE], E, c, h = (x - x0)/1000000, global_err = 0.;
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

int euler(double x0, double x, double* y0, double* y, int size, int num_steps) {
	int i, j;
	double h = (x - x0)/num_steps, buff[ST_SIZE];
	memcpy(y, y0, size * sizeof(double));
	for(i = 0; i < num_steps; i++, x0+=h) {
		f(x0, y, size, buff);
		for(j = 0; j < size; j++)
			y[j] += h*buff[j];
	}
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

void plot(const char* fname) { 
        FILE * gnuplot_pipe = popen("gnuplot -persistent", "w");
        fprintf(gnuplot_pipe, "%s '%s' %s\n", "plot", fname, "with line");
}


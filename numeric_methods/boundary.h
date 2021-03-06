#pragma once

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

double** create_matrix(int rows, int cols);
void delete_matrix(double** A, int rows);
int read_matrix(double** m, int rows, int cols, FILE* inp);
void print(double** m, int rows, int cols, FILE* out);
void print_vector(double* y, int size, FILE* out);
double norm(double* y, int size, const char name);
int multiply(double** A, double** B, double** res, int arows, int acols, int brows, int bcols);
void triangle(double** m, double* b, int n);
void gauss(double** m, double* x, double* b, int n);

double chord_method(double x1, double x2, double* y1, double* y2, int size);
int euler(double x0, double x, double* y0, double* y, int size, int num_steps);
double dormand5(double x0, double* y, double* y1, int size, double h);
double butcher(double x0, double* y, double* y1, int size, double h);
double ingland(double x0, double* y, double* y1, int size, double h);
double dormand8(double x0, double* y, double* y1, int size, double h);
double integrate(double x0, double x, double* y0, double *y, int size, int num_steps, 
		double iteration(double x0, double* y, double* y1, int size, double h));
double integrate_autostep(double x0, double x, double* y0, double* y, int size, double err_min, double err_max, double h, 
		double iteration(double x0, double* y, double* y1, int size, double h), int flag);

void runge_numbers(double x0, double x, double* y0, int size);
double track(double a, double b, double* y0, int size, int num_points, int plot);

int shoot(double a, double b, double* y0, int size, int k, double eps, void discrepancy(double* y0, double* y, double* res), double err_min, double err_max);
double inception(double* p, double alpha1, double alpha2, int num_steps, double a, double b, double* y0, int size, int k, double eps,
	       	void discrepancy(double* y0, double* y, double* res), double err_min, double err_max);


double eigen_value(double x, double* y);
void discrepancy(double* y0, double* y, double* res);
void start_value(double* y0);

void derivative(int n, double* x, void f(double* x, double* y), double** res, double* buff);
int newton_method(int n, double* x, void f(double* x, double* y), double eps);

#define PI 3.141592653589793238462643
#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b)
#define abs(a) a > 0 ? a : -1*a
#define random(a, b) ((double)rand())*((b) - (a))/(double)RAND_MAX + (a)
#define is_zero(a) ((a) < 1.e-55 && (a) > -1.e-55)
#define DELTA 1.e-8
#define ST_SIZE 8

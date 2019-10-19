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

int euler(double x0, double x, double* y0, double* y, int size, int num_steps);
double runge(double x0, double x, double *y0, double *y, int size, int num_steps);
double runge_with_autostep(double x0, double x, double* y0, double* y, int size, double err_min, double err_max);
double runge_hardcore(double x0, double x, double* y0, double* y, int size, double err_min, double err_max);
void runge_numbers(double x0, double x, double* y0, int size);
void plot(const char* fname);

int shooting(double a, double b, double* y0, int size, int k, double eps);

void f(double x, double* y, int size, double* res);
double eigen_value(double x, double* y);
void discrepancy(double* y, double* v);
void jacobian(double** m);

#define MAX(a,b) (a)>(b)?(a):(b)
#define is_zero(a) ((a) < 0.000000001 && (a) > -0.000000001)
#define DELTA 1.e-8
#define ST_SIZE 8

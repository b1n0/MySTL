#pragma once

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

double** create_matrix(int rows, int cols);
void delete_matrix(double** A, int rows);
int read_matrix(double** m, int rows, int cols, FILE* inp);
void print(double** m, int rows, int cols);
void print_vector(double* y, int size, FILE* out);
double norm(double* y, int size);
int multiply(double** A, double** B, double** res, int arows, int acols, int brows, int bcols);


void f(double x, double* y, int size, double* res);
int euler(double x0, double x, double* y0, double* y, int size, int num_steps);
int runge(double x0, double x, double *y0, double *y, int size, int num_steps);
int runge_with_autostep(double x0, double x, double* y0, double* y, int size, double err_min, double err_max);
int gauss(double** m, double* x, double* b, int n);
int shooting(double* y0, int size, int k, double a, double b);

#define ST_SIZE 32

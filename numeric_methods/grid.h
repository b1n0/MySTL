
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
void running_method(int n, double* x, double** c, double* f);

double** create_matrix(int rows, int cols);
void delete_matrix(double** A, int rows);
void explicit_solve(double** u, int tn, int xn, double a);
void implicit_solve(double** u, int tn, int xn, double a);
void inaccuracy(double* u, double* v, int xn, double a, double* ans);
void save(double** u, int tn, int xn, const char* fname);
void plot(double** u, int tn, int xn, const char* fname);

#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b)
#define abs(a) a > 0 ? a : -1*a
#define is_zero(a) ((a) < 1.e-55 && (a) > -1.e-55)

#define Tn 10001
#define Xn 10001

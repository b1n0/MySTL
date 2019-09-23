#pragma once

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

void f(double x, double* y, int size, double* res);
int euler(double x0, double x, double* y0, double* y, int size, int num_steps);
int runge(double x0, double x, double *y0, double *y, int size, int num_steps);
int runge_with_autostep(double x0, double x, double* y0, double* y, int size, double h, double err, double K);
void print(double* y, int size);
int plot(double* y, int num_points);
int gauss(double** m, double* x, double* b, int n);

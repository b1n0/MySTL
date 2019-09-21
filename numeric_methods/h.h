#pragma once

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

void f(double x, double* y, int size, double* res);
int runge(double x0, double x, double *y0, double *y, int size, int num_steps);
int runge_with_autostep(double x0, double x, double* y0, double* y, int size, double h, double err, double K);
int euler(double x0, double x, double* y0, double* y, int size, int num_steps);
void print(double* y, int size);
int plot(double* y, int num_points);
int derivative(double x, double* y, double* res, int size, double h);

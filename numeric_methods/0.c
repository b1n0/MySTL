<<<<<<< HEAD
#include "h.h"

int derivative(double x, double* y, double* res, int size, double h) {
	int i;
	double buff[32];
	f(x, y, size, res);
	f(x + h, y, size, buff);
	for(i = 0; i < size; i++)
		res[i] = (buff[i] - res[i])/h;
	return 0;
}
void f(double x, double* y, int size, double* res) {
       res[0] = (exp(2.*x) - 1)/(exp(2.*x) + 1);
}       

int main(void) {
	int i;
	double y[1], res[1];
	derivative(0, y, res, 1, 0.01);	
	print(res,1);
	return 0;
}


=======
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
double f(double x);
double f(double x) { return (exp(2*x) - 1)/(exp(2*x) + 1); } 
	
int main(void) {
	double x = 0, h = 1;
	for (int i = 1; i < 30; i++) {
		h *= 0.1;
		printf("%25.20lf\n", ((f(x+h) - f(x))/h) - 1.);
	}
	return 0;
}
>>>>>>> a68f61229313779aff5f5379be5325224cd15ec0

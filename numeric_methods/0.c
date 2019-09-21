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

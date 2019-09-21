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



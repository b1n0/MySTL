#include "h.h"

void f(double x, double* y, int size, double* res) {
	int i;
	for(i = 0; i < size; i++) res[i] = x;
}

int main(void) {
	double x0 = 0, x = 1, y0 = 0, y = 0;
	evaluate(x0, x, &y0, &y, 1, 100);	
	printf("%lf\n", y);
	return 0;
}

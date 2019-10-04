#include <stdlib.h>
#include <stdio.h>
#include <math.h>
double up(int k);
double down(int k);
double integrate(int k, int num_steps);

int main(void) {
	printf("%lf\n", up(30));
	printf("%lf\n", down(30));
	printf("%lf\n", integrate(30, 1000000));
	return 0;
}

double up(int k) {
	double res = log(7./6.);
	for(int n = 1; n <= k; n++)
		res = 1./n - 6.*res;
	return res;
}

double down(int k) {
	double res = 0;
	for(int n = 100; n >= k; n--)
		res = 1./(6.*n) - res/6.;
	return res;
}

double integrate(int k, int num_steps) {
	double res = 0, h = 1./num_steps, x = 0;
	for(int i = 0; i < num_steps; i++, x+=h)
		res += h*(pow(x, k)/(x + 6));
	return res;
}

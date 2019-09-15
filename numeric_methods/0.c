#include "h.h"

int derivative(x, y, res, size, h) {
	int i;
	double buff[ST_SIZE];
	f(x, y, size, res);
	f(x + h, buff);
	for(i = 0; i < size; i++)
		res[i] = (buff[i] - res[i])/h;
	return 0;
}



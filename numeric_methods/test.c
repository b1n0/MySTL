#include "h.h"

int main(void) {
	double **m = create_matrix(1, 1);
	double **x = create_matrix(1, 1);
	double **b = create_matrix(1, 1);
	read_matrix(m, 1, 1, stdin);
	read_matrix(b, 1, 1, stdin);
	printf("\n");
	
	gauss(m, x[0], b[0], 1);
	print(m, 1, 1);
	print_vector(x[0], 1, stdout);

	delete_matrix(m, 1);
	delete_matrix(x, 1);
	delete_matrix(b, 1);
	return 0;
}

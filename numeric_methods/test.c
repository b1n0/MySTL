#include "h.h"

int main(void) {
	double **m = create_matrix(3, 3);
	double **x = create_matrix(1, 3);
	double **b = create_matrix(1, 3);
	read_matrix(m, 3, 3, stdin);
	read_matrix(b, 1, 3, stdin);
	printf("\n");
	
	gauss(m, x[0], b[0], 3);
	multiply(x, m, b, 3, 1, 3, 3);
	print_vector(x[0], 3, stdout);
	printf("\n");
	print_vector(b[0], 3, stdout);

	delete_matrix(m, 3);
	delete_matrix(x, 3);
	delete_matrix(b, 3);
	return 0;
}

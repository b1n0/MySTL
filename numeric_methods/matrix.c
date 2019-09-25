#include "h.h"
#include <string.h>

double** create_matrix(int rows, int cols) {
	double** A;
	A = (double**)malloc(rows*sizeof(double*));
	for(int i = 0; i < rows; i++)
		A[i] = (double*)malloc(cols*sizeof(double));
	return A;
}

void delete_matrix(double** A, int rows) {
	for(int i = 0; i < rows; i++)
		free(A[i]);
	free(A);
}

int read_matrix(double** m, int rows, int cols, FILE* inp) {
	if (m) {
		int i,j,l;
		char str[2048], *s;
		for(i = 0; i < rows; i++) {
			fgets(str, 2048, inp);
			for(s = str, l = 0, j = 0; j < cols; j++, s+=l)
				sscanf(s, "%lf%n", m[i]+j, &l);
		}
		return 0;
	}
	return -1;
}
void print(double** m, int rows, int cols) {
	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++)
			printf("%lf ", m[i][j]);
		printf("\n");
	}
}
int multiply(double** A, double** B, double** res, int arows, int acols, int brows, int bcols) {
	if(A && B && res && acols == brows) {
		for(int i = 0; i < arows; i++) {
			for(int j = 0; j < bcols; j++) {
				res[i][j] = 0.;
				for(int k = 0; k < brows; k++)
					res[i][j] += A[i][k] * B[k][j];
			}
		}
		return 0;
	}
	return -1;
}

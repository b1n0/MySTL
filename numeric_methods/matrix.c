#include "h.h"
#define is_zero(a) ((a) < 0.000000001 && (a) > -0.000000001)

int gauss(double** m, double* x, double* b, int n) {
        int i, j, k, ind, c = 0, rank = 0;
        double* tmp;
        for(k = 0; k < n; k++) {
                for(ind = -1, i = c; i < n; i++)
                        if (!is_zero(m[i][k])) {
                                ind = i;
                                for(j = 0; j < k; j++)
                                        m[i][j] /= m[i][k];
                                b[i] /= m[i][k];
                        }
                for(i = c; i < n && ind != -1; i++) {
                        if (!is_zero(m[i][k]) && i != ind) {
                                for(j = k; j < n; j++)
                                        m[i][j] -= m[ind][j];
                                b[i] -= b[ind];
                        }
                }
                if (ind != -1 && ind != c) {
                        tmp = m[ind]; m[ind] = m[c]; m[c] = tmp;
                        b[ind] += b[c]; b[c] = b[ind] - b[c]; b[ind] = b[ind] - b[c];
                }
                c++;
        }

        for(k = 0; k < n; k++) {
                double s;
                if is_zero(m[n-1-k][n-1-k])  x[n-1-k] = 0;
                else {
                        for(s = 0., i = 0; i < k; i++)
                                s-= x[n-1-i]*m[n-1-k][n-1-i];
                        x[n-1-k] = (b[n-1-k] + s)/m[n-1-k][n-1-k];
                }
        }
        return 0;
}

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




#include "h.h"

void gauss(double** m, double* x, double* b, int n) {
	double s = 0.;
	triangle(m, b, n);
        for(int k = 0; k < n; k++) {
                if is_zero(m[n-1-k][n-1-k])  
			x[n-1-k] = 0;
                else {
                        for(int i = 0, s = 0.; i < k; i++)
                                s-= x[n-1-i]*m[n-1-k][n-1-i];
                        x[n-1-k] = (b[n-1-k] + s)/m[n-1-k][n-1-k];
                }
        }
}

void triangle(double** m, double* b, int n) {
	int i, j, k, ind_nz, c;
	double* tmp;
	for(ind_nz = -1, c = 0, k = 0; k < n; k++, c++) {	
		for(ind_nz = -1, i = c; i < n; i++) {
			if(!is_zero(m[i][k])) {
				ind_nz = ind_nz == -1 ? i : ind_nz;
				b[i] /= m[i][k];	
				for(j = n-1; j >= k; j--)
				       m[i][j] /= m[i][k];
			}
		}
		if(ind_nz != -1) {
			for(i = c; i < n; i++) {
				if(!is_zero(m[i][k]) && ind_nz != i) {
					for(j = k; j < n; j++)
						m[i][j] -= m[ind_nz][j];
					b[i] -= b[ind_nz];	
				}
			}	
			if(ind_nz != c) {
                        	tmp = m[ind_nz]; m[ind_nz] = m[c]; m[c] = tmp;
                        	b[ind_nz] += b[c]; b[c] = b[ind_nz] - b[c]; b[ind_nz] = b[ind_nz] - b[c];
			}
		}	
	}
}

double** create_matrix(int rows, int cols) {
	double** A;
	A = (double**)malloc(rows*sizeof(double*));
	for(int i = 0; i < rows; i++) {
		A[i] = (double*)malloc(cols*sizeof(double));
		for(int j = 0; j < cols; j++) A[i][j] = 0.;
	}
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

void print(double** m, int rows, int cols, FILE* out) {
	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++)
			fprintf(out, "%lf ", m[i][j]);
		printf("\n");
	}
}

void print_vector(double* y, int size, FILE* out) {
	for(int i = 0; i < size; i++)
		fprintf(out, "%lf ", y[i]);
	fprintf(out, "\n");
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

double norm(double* y, int size, const char name) {
	double l = 0.;
	if(name == 'e') {
		for(int i = 0; i < size; i++)
			l += y[i]*y[i];
		l = sqrt(l);
	}
	else if (name == 'm') 
		for(int i = 0; i < size; i++)
			l += fabs(y[i]);
	return l;
}

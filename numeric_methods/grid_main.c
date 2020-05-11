#include "h.h"

int main(void) {
	int i, j, k, xn, tn;
	double ans[4], **u, **v, a;
	v = create_matrix(Tn, Xn);
	for(k = 0, a = 1.; k < 2; k++, a = 0., printf("\n")) {
		implicit_solve(v, Tn, Xn, a);
		plot(v, Tn, Xn, "plot.txt");
		for(i = 0, tn = 11; i < 3; i++, tn = (tn-1)*10 + 1)
			for(j = 0, xn = 11; j < 3; j++, xn = (xn-1)*10 + 1) { 
				u = create_matrix(tn, xn);
				implicit_solve(u, tn, xn, a);
				inaccuracy(u[tn-1], v[Tn-1], xn, a, ans);
				printf("%lf| %lf| %lf %lf %lf %lf \n",1./(tn-1), 1./(xn-1), ans[0], ans[1], ans[2], ans[3]);
				delete_matrix(u, tn);
			}	
	}
	delete_matrix(v, Tn);
	return 0;
}


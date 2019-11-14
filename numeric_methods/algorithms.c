#include "h.h"

double chord_method(double x1, double x2, double* y1, double* y2, int size) {
	double y[ST_SIZE], xroot;
	if(g(y1)*g(y2) < 0) {
		while(fabs(g(y1)) > 1.e-15 && fabs(g(y2)) > 1.e-15) {
			xroot = x1 - g(y1)*(x2 - x1)/(g(y2) - g(y1));
			dormand8(x1, y1, y, size, xroot - x1);
			if(g(y1)*g(y) > 0) { memcpy(y1, y, size*sizeof(double)); x1 = xroot; }
			else { memcpy(y2, y, size*sizeof(double)); x2 = xroot; }
		}
		return fabs(g(y1)) > fabs(g(y2)) ? x2 : x1;
	}
	return x2;
}

double integrate(double x0, double x, double* y0, double *y, int size, int num_steps,
		double iteration(double x0, double* y, double* y1, int size, double h)) {
        double h = (x - x0)/num_steps, err, global_err = 0.;
	memcpy(y, y0, size*sizeof(double));
        for (int i = 0; i < num_steps; i++, x0+=h) {
		err = iteration(x0, y, y, size, h);
		global_err = exp(h*eigen_value(x, y))*global_err + err;
        }
        return global_err;
}

double integrate_autostep(double x0, double x, double* y0, double *y, int size, double err_min, double err_max, double h, 
		double iteration(double x0, double* y, double* y1, int size, double h), int flag) {
        int i; 
	double y1[ST_SIZE], h_old, xroot, global_err = 0., err = err_max + 1., x_max = 0., y_max = 0.;
        double fac = 0.8, facmin = 0.2, facmax = 2.;
	memcpy(y, y0, size*sizeof(double));
	for(i = 0; (x0 < x - h && h > 0) || (x0 > x - h && h < 0); i++) {
	//	while(err >= err_max) {
		while(1) { 
			err = iteration(x0, y, y1, size, h);
			h_old = h;
			if(!is_zero(err)) h *= MIN(facmax, MAX(facmin, fac*pow(err_max/err, 1./8.)));	
			else h *= err < err_min ? 2. : 1.;
	//	}
			if(err < err_max) {
				x0 = chord_method(x0, x0 + h_old, y, y1, size);
				break;
			}
		}
	//	x0 = chord_method(x0, x0 + h_old, y, y1, size);
		global_err = exp(h_old*eigen_value(x, y1))*global_err + err; 
		memcpy(y, y1, size*sizeof(double));
		x_max = MAX(x_max, fabs(y[0] - sin(x0+h_old)));
		y_max = MAX(y_max, fabs(y[1] - cos(x0+h_old)));
	}
	if(flag) printf("%d | %1.15lf | %1.15lf | ", i, x_max, y_max);
	return global_err + integrate(x0, x, y, y, size, 10, iteration);
}

double dormand5(double x0, double* y, double* y1, int size, double  h) {
	double k[7][ST_SIZE], buff[ST_SIZE], err = 0.;
	int i;
	f(x0, y, k[0]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*k[0][i]*0.2;
	f(x0 + 0.2*h, buff, k[1]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(3./40.) + k[1][i]*(9./40.));
	f(x0 + 0.3*h, buff, k[2]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(44./45.) - k[1][i]*(56./15.) + k[2][i]*(32./9.));
	f(x0 + 0.8*h, buff, k[3]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(19372./6561.) - k[1][i]*(25360./2187.) + k[2][i]*(64448./6561.) - k[3][i]*(212./729.));
	f(x0 + (8./9.)*h, buff, k[4]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(9017./3168.) - k[1][i]*(355./33.) + k[2][i]*(46732./5247.) + k[3][i]*(49./176.) - k[4][i]*(5103./18656.));
	f(x0 + h, buff, k[5]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(35./384.) + k[2][i]*(500./1113.) + k[3][i]*(125./192.) - k[4][i]*(2187./6784.) + k[5][i]*(11./84.));
	f(x0 + h, buff, k[6]);

	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(5179./57600.) + k[2][i]*(7571./16695.) + k[3][i]*(393./640.) - k[4][i]*(92097./339200.) + k[5][i]*(187./2100.) + k[6][i]*(1./40.));
	for(i = 0; i < size; i++) y1[i] = y[i] + h*(k[0][i]*(35./384.) + k[2][i]*(500./1113.) + k[3][i]*(125./192.) - k[4][i]*(2187./6784.) + k[5][i]*(11./84.));

	for(i = 0; i < size; i++) err += fabs(y1[i] - buff[i]);
       	return err;
}

double butcher(double x0, double* y, double* y1, int size, double  h) {
	double k[7][ST_SIZE], buff[ST_SIZE], y2[ST_SIZE], err = 0.;
	int i;
	f(x0, y, k[0]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*k[0][i]*0.5;
	f(x0 + h*0.5, buff, k[1]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(2./9.) + k[1][i]*(4./9.));
	f(x0 + h*(2./3.), buff, k[2]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(7./36.) + k[1][i]*(2./9.) - k[2][i]*(1./12.));
	f(x0 + h*(1./3.), buff, k[3]);
	for(i = 0; i < size; i++) buff[i] = y[i] - h*(k[0][i]*(35./144.) - k[1][i]*(55./36.) + k[2][i]*(35./48.) + k[3][i]*(15./8.));
	f(x0 + h*(5./6.), buff, k[4]);
	for(i = 0; i < size; i++) buff[i] = y[i] - h*(k[0][i]*(1./360.) - k[1][i]*(11./36.) - k[2][i]*0.125 + k[3][i]*0.5 + k[4][i]*0.1);
	f(x0 + h*(1./6.), buff, k[5]);
	for(i=0;i<size;i++) buff[i]=y[i] + h*(-k[0][i]*(41./260.) + k[1][i]*(22./13.) + k[2][i]*(43./156.) - k[3][i]*(118./39.) + k[4][i]*(32./195.)+k[5][i]*(80./39.));
	f(x0 + h, buff, k[6]);

	for(i = 0; i < size; i++) y2[i] = y[i]+(k[0][i]*13./200. + k[2][i]*11./40. + k[3][i]*11./40. + k[4][i]*0.16 + k[5][i]*0.16 + k[6][i]*13./200.)*h;
	
	f(x0, y, k[0]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*k[0][i]*0.5;
	f(x0 + h*0.5, buff, k[1]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i] + k[1][i])*0.25;
	f(x0 + h*0.5, buff, k[2]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(2*k[2][i] - k[1][i]);
	f(x0 + h, buff, k[3]);	
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(7*k[0][i] +10*k[1][i] + k[3][i])/27;
	f(x0 + h*(2./3.), buff, k[4]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(28*k[0][i] - 125*k[1][i] + 546*k[2][i] + 54*k[3][i] - 378*k[4][i])/625;
	f(x0 + h*0.2, buff, k[5]);
	
	for(i = 0; i < size; i++) y1[i] = y[i] + (14*k[0][i] + 35*k[3][i] + 162*k[4][i] + 125*k[5][i])*h/336.;

	for(i = 0; i < size; i++) err += fabs(y1[i] - y2[i]);
	return err;
}

double ingland(double x0, double* y, double* y1, int size, double  h) {
	double k[6][ST_SIZE], buff[ST_SIZE], err = 0.;
	int i;
	f(x0, y, k[0]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*k[0][i]*0.5;
	f(x0 + h*0.5, buff, k[1]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i] + k[1][i])*0.25;
	f(x0 + h*0.5, buff, k[2]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(2.*k[2][i] - k[1][i]);
	f(x0 + h, buff, k[3]);	
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(7.*k[0][i] +10.*k[1][i] + k[3][i])*(1./27.);
	f(x0 + h*(2./3.), buff, k[4]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(28.*k[0][i] - 125.*k[1][i] + 546.*k[2][i] + 54.*k[3][i] - 378.*k[4][i])*(1./625.);
	f(x0 + h*0.2, buff, k[5]);

	for(i = 0; i < size; i++) y1[i] = y[i] + (14.*k[0][i] + 35.*k[3][i] + 162.*k[4][i] + 125.*k[5][i])*h/336.;
	for(i = 0; i < size; i++)  err += fabs(((-42.)*k[0][i] - 244.*k[2][i] - 21.*k[3][i] + 162.*k[4][i] + 125.*k[5][i])*h*(1./336.));

	return err;
}

double dormand8(double x0, double* y, double* y1, int size, double  h) {
	double k[13][ST_SIZE], buff[ST_SIZE], err = 0.;
	int i;	
	f(x0, y, k[0]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*k[0][i]*(1./18.);
	f(x0 + h*(1./18.), buff, k[1]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(1./48.) + k[1][i]*(1./16.));
	f(x0 + h*(1./12.), buff, k[2]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(1./32.) + k[2][i]*(3./32.));
	f(x0 + h*0.125, buff, k[3]);	
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(5./16.) - k[2][i]*(75./64.) + k[3][i]*(75./64.));
	f(x0 + h*(5./16.), buff, k[4]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(3./80.) + k[3][i]*(3./16.) + k[4][i]*(3./20.));
	f(x0 + h*(3./8.), buff, k[5]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(29443841./614563906.) + k[3][i]*(77736538./692538347.) - k[4][i]*(28693883./1125000000.) + k[5][i]*(23124283./1800000000.));
	f(x0 + h*(59./400.), buff, k[6]);
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(16016141./946692911.) +  k[3][i]*(61564180./158732637.) + k[4][i]*(22789713./633445777.) + k[5][i]*(545815736./2771057229.) - k[6][i]*(180193667./1043307555.));
	f(x0 + h*(93./200.), buff, k[7]); 
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(39632708./573591083.) - k[3][i]*(433636366./683701615.) - k[4][i]*(421739975./2616292301.) + k[5][i]*(100302831./723423059.) + k[6][i]*(790204164./839813087.) + k[7][i]*(800635310./3783071287.));
	f(x0 + h*(5490023248./9719169821.), buff, k[8]); 
	for(i = 0; i < size; i++) buff[i] = y[i] + h*(k[0][i]*(246121993./1340847787.) - k[3][i]*(37695042795./15268766246.) - k[4][i]*(309121744./1061227803.) - k[5][i]*(12992083./490766935.) + k[6][i]*(6005943493./2108947869.) + k[7][i]*(393006217./1396673457.) + k[8][i]*(123872331./1001029789.));
	f(x0 + h*(13./20.), buff, k[9]); 
	for(i = 0; i < size; i++) buff[i] = y[i]+h*(-k[0][i]*(1028468189./846180014.) + k[3][i]*(8478235783./508512852.) + k[4][i]*(1311729495./1432422823.) - k[5][i]*(10304129995./1701304382.) - k[6][i]*(48777925059./3047939560.) + k[7][i]*(15336726248./1032824649) - k[8][i]*(45442868181./3398467696.) + k[9][i]*(3065993473./597172653.));
	f(x0 + h*(1201146811./1299019798.), buff, k[10]); 
	for(i = 0; i < size; i++) buff[i] = y[i]+h*(k[0][i]*(185892177./718116043.) - k[3][i]*(3185094517./667107341.) - k[4][i]*(477755414./1098053517.) - k[5][i]*(703635378./230739211.) + k[6][i]*(5731566787./1027545527.) + k[7][i]*(5232866602./850066563.) - k[8][i]*(4093664535./808688257.) + k[9][i]*(3962137247./1805957418.) + k[10][i]*(65686358./487910083.));
	f(x0 + h, buff, k[11]); 
	for(i = 0; i < size; i++) buff[i] = y[i]+h*(k[0][i]*(403863854./491063109.) - k[3][i]*(5068492393./434740067.) - k[4][i]*(411421997./543043805.) + k[5][i]*(652783627./914296604.) + k[6][i]*(11173962825./925320556.) - k[7][i]*(13158990841./6184727034.) + k[8][i]*(3936647629./1978049680.) - k[9][i]*(160528059./685178525.) + k[10][i]*(248638103./1413531060.));
	f(x0 + h, buff, k[12]); 

	for(i = 0; i < size; i++) buff[i] = y[i] + (k[0][i]*(14005451./335480064.) - k[5][i]*(59238493./1068277825.) + k[6][i]*(181606767./758867731.) + k[7][i]*(561292985./797845732.) - k[8][i]*(1041891430./1371343529.) + k[9][i]*(760417239./1151165299.) + k[10][i]*(118820643./751138087.) - k[11][i]*(528747749./2220607170.) + k[12][i]*0.25)*h;
	for(i = 0; i < size; i++) y1[i] = y[i] + (k[0][i]*(13451932./455176623.) - k[5][i]*(808719846./976000145.) + k[6][i]*(1757004468./5645159321.) + k[7][i]*(656045339./265891186.)- k[8][i]*(3867574721./1518517206.) + k[9][i]*(465885868./322736535.) +	k[10][i]*(53011238./667516719.) + k[11][i]*(2./45.))*h;

	for(i = 0; i < size; i++) err += fabs(y1[i] - buff[i]);
	//for(i = 0; i < size; i++) err += fabs(h*(0.012194277465176744*k[0][i] + 0.7731539478765578*k[5][i] - 0.07192809284993823*k[6][i] - 1.763834521196444*k[7][i] + 1.787182038027448*k[8][i] - 0.7829855527544889*k[9][i] + 0.07877188662899604*k[10][i] + 0.19366509430841836*k[11][i] + 0.25*k[12][i]));
       	return err;
}

int euler(double x0, double x, double* y0, double* y, int size, int num_steps) {
	int i, j;
	double h = (x - x0)/num_steps, buff[ST_SIZE];
	memcpy(y, y0, size * sizeof(double));
	for(i = 0; i < num_steps; i++, x0+=h) {
		f(x0, y, buff);
		for(j = 0; j < size; j++)
			y[j] += h*buff[j];
	}
	return 0;
}

void runge_numbers(double x0, double x, double* y0, int size) {
        double y8[ST_SIZE], y10[ST_SIZE], y12[ST_SIZE];
        integrate_autostep(x0, x, y0, y8, size, 1.e-9, 1.e-8, 1.e-2, dormand8, 0);
        integrate_autostep(x0, x, y0, y10, size, 1.e-11, 1.e-10, 1.e-2, dormand8, 0);
        integrate_autostep(x0, x, y0, y12, size, 1.e-13, 1.e-12, 1.e-2, dormand8, 0);
        printf("runge number in %lf \t ", x);
        for(int i = 0; i < size; i++) printf("%lf ", (y8[i] - y10[i])/(y10[i] - y12[i]));
        printf("\n");
}

double track(double a, double b, double* y0, int size, int num_points, int plot) { 
	int i,j;
	FILE *gnuplot_pipe, *f = fopen("plt.txt", "w");
	double y[ST_SIZE], x, h = (b - a)/num_points, global_err = 0.;
	memcpy(y, y0, size*sizeof(double));
	for(x = a, i = 0; i < num_points; i++, x+=h) {
		for(j = 0, fprintf(f, "%lf\t", x); j < size; j++) fprintf(f, "%lf\t", y[j]);
		global_err += integrate_autostep(x, x+h, y, y, size, 1.e-9, 1.e-8, 1.e-6, dormand8, 0);
		fprintf(f, "%lf\n", u(y));
	}
	for(j = 0, fprintf(f, "%lf\t", x); j < size; j++) fprintf(f, "%lf\t", y[j]);
	fprintf(f, "%lf", u(y));
	fclose(f);
	if(plot) {
		gnuplot_pipe = popen("gnuplot -persistent", "w");
        	fprintf(gnuplot_pipe, "%s%d%s\n", "set key off; plot for[col=2:", size + 1, ":1]'plt.txt' using 1:col with lines");
	}
	return global_err;
}


#define NUM_F_MODES 9

void print_fourier(unsigned long long step)
{
int ii,kk;
double delta_r_sh[DIMENSION];//(R-Rsh)
double theta, phi;
double rad, xyrad;
double inv_n=1./num_shell_monos;

vector<double> re_fourier_tot(NUM_F_MODES, 0.);
vector<double> im_fourier_tot(NUM_F_MODES, 0.);

for(ii = NUMBER_IN_POLYMER; ii < NUMBER_IN_POLYMER + num_shell_monos; ii++)
{
	rad=0.;
	for(kk = 0; kk < DIMENSION; kk++)
	{
	  delta_r_sh[kk] = mono_list[ii].get_prev_pos(kk) - shell_cm[kk];
	  rad += delta_r_sh[kk]*delta_r_sh[kk];
	}
	xyrad = rad-delta_r_sh[0]*delta_r_sh[0];
	rad = sqrt(rad);
	xyrad = sqrt(xyrad);
	

#if PRINT_FOURIER
	for(kk = 0; kk < NUM_F_MODES; kk++)
	{
	 re_fourier_tot[kk] += (xyrad*re_exp_nphi(delta_r_sh[1], delta_r_sh[2], xyrad, kk));
	 im_fourier_tot[kk] += (xyrad*im_exp_nphi(delta_r_sh[1], delta_r_sh[2], xyrad, kk));
	}

#endif//fourier
}//for(ii loop over shell monos



#if PRINT_FOURIER
fprintf(fourierfile, "%llu ", step);
for(kk = 0; kk < NUM_F_MODES; kk++)
{
        re_fourier_tot[kk] *= (INV_PI*inv_n);
        im_fourier_tot[kk] *= (INV_PI*inv_n);

	fprintf(fourierfile, "%g %g ", re_fourier_tot[kk], im_fourier_tot[kk]);
}

fprintf(fourierfile, "\n");
fflush(fourierfile);
#endif

}//end print_fourier()



double re_exp_nphi(double x, double y, double r, int n)
{
if(n > 1)
{
if(x>0.)
	return cos(((double)n)*atan(y/x));
else
	return cos(((double)n)* (atan(y/x)+PI));
}
else if(n == 1)
{
	return x/r;
}
else//n == 0
{
	return 1.;
}
}


double im_exp_nphi(double x, double y, double r, int n)
{
if(n > 1)
{
if(x>0.)
        return sin(((double)n)*atan(y/x));
else
        return sin(((double)n)*(atan(y/x)+PI));
}
else if(n == 1)
{
        return y/r;
}
else//n == 0
{
        return 0.;
}
}







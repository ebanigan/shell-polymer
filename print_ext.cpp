void print_ext(unsigned long long step)
{
     char filename[96];
     FILE *outfile;
     long ii;
     int kk;
     double total[DIMENSION];
     double temp_pos[DIMENSION];
     int dummy;     

if(step >= 0)
{

/********calculate CM and print if desired*********/
for(ii = 0; ii < NUMBER_IN_POLYMER; ii++)
{
 for(kk = 0; kk < DIMENSION; kk++)
 {
  temp_pos[kk] = mono_list[ii].get_prev_pos(kk);
  total[kk] += temp_pos[kk];
 }
}

if(NUMBER_IN_POLYMER > 0)
{
for(kk = 0; kk < DIMENSION; kk++)
  center_mass_pos[kk] = total[kk]/NUMBER_IN_POLYMER;

#if PRINT_CM_ONLY            
     sprintf(filename, "output/polymer_move%6.6i", TRIALNUMBER);
    
     if(step == 0)
      outfile = fopen(filename, "w");
     else
      outfile = fopen(filename, "a");
#endif
}
/********************************************/
}//if step >=0




/*shell rgyr (in x direction)*/
double shell_rgx2 = 0.;
double min_shell = 9.*L[0];
double max_shell = -9.*L[0];

double max_y = -9.*L[1];
double max_z = -9.*L[2];
double min_y = 9.*L[1];
double min_z = 9.*L[2];
double monox, monoy, monoz;

double tot_spring_extension = 0.;

//SHELL RG
for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
	monox = mono_list[ii].get_prev_pos(0);
	monoy = mono_list[ii].get_prev_pos(1);
	monoz = mono_list[ii].get_prev_pos(2);


	if(monox < min_shell)
		min_shell = monox; 
	else if(monox > max_shell)
		max_shell = monox;

	if(monoy < min_y)
		min_y = monoy;
	else if(monoy > max_y)
		max_y = monoy;

	if(monoz < min_z)
		min_z = monoz;
	else if(monoz > max_z)
		max_z = monoz;

	shell_rgx2 = shell_rgx2 + (monox - shell_cm[0])*(monox - shell_cm[0]);

	if(LENGTH_CONTROLLED_LOAD)
        if(mono_list[ii].get_load_mono())
	{
	   tot_spring_extension += mono_list[ii].get_sign_load()*(0.5*L[0] + mono_list[ii].get_sign_load()*(pipette_dist_from_center - mono_list[ii].get_load_rest_length()) - mono_list[ii].get_prev_pos(0));
	}

}

if(step >= 0)
{
if(LENGTH_CONTROLLED_LOAD)
{
  fprintf(totforcefile, "%llu %g %g\n", step, PIPETTE_STIFFNESS*tot_spring_extension, pipette_dist_from_center);
  fflush(totforcefile);
}

if(num_shell_monos > 0)
{
if(F_LOAD < -1.e-12)
fprintf(extfile, "%llu %g %g %g %g %g\n", step, -2., max_y-min_y, max_z-min_z, shell_rgx2/num_shell_monos, mono_list[indentation2].get_prev_pos(0) - mono_list[indentation1].get_prev_pos(0));
else
fprintf(extfile, "%llu %g %g %g %g %g\n", step, -2., max_y-min_y, max_z-min_z, shell_rgx2/num_shell_monos, max_shell-min_shell);
}
else
{
fprintf(extfile, "%llu %g %g %g %g %g\n", step, 0., 0., 0., 0., 0.);
}

fflush(extfile);
}//if step>=0

//movie//if(false)
if(BOX_RESIZE)
if(step % (NUMSKIP*100) == 0)
{

 resize_box(max_shell-min_shell, max_y-min_y, max_z-min_z);
  

 if(!THERMAL)
 {
  if((step > 15500000) && (fabs(max_shell-min_shell - saved_extension) < (0.05*SHELL_RADIUS)*1.e-6)) 
  {
	fprintf(stderr, "converged, exiting\n");

        char mvcommand[128];
        sprintf(mvcommand, "mv /tmp/[rs]*%6.6i output/", TRIALNUMBER);
        dummy= system(mvcommand);

	exit(1);
  }
  else
  {
   saved_extension = max_shell-min_shell;
  }

 }
 
 //saved_extension = 0.;
}
/*
else 
{
	saved_extension += (max_shell-min_shell)*2./99.;
}
*/

#if PRINT_CM_ONLY
if(NUMBER_IN_POLYMER > 0)
{
fprintf(outfile, "%li %g %g %g\n", step, center_mass_pos[0], center_mass_pos[1], center_mass_pos[2]);

     fflush(outfile);
     fclose(outfile);
}
#endif



}


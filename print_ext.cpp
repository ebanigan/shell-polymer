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
double monox;

double tot_shell_force=0.;
vector<double> extensions;
double shell_extension=0.;
double yspan, zspan;

//SHELL RG
for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
        monox = mono_list[ii].get_prev_pos(0);
	shell_rgx2 = shell_rgx2 + (monox - shell_cm[0])*(monox - shell_cm[0]);
}

//force
if(LENGTH_CONTROLLED_LOAD) // calculate instantaneous force on pipette / nucleus
  tot_shell_force=compute_applied_force();

extensions = compute_extension();
shell_extension= extensions[0];
yspan=extensions[1];
zspan=extensions[2];

if(step >= 0)
{
if(LENGTH_CONTROLLED_LOAD)
{
  fprintf(totforcefile, "%llu %g %g\n", step, tot_shell_force, pipette_dist_from_center);
  fflush(totforcefile);
}

if(num_shell_monos > 0)
{
if(F_LOAD < -1.e-12)
fprintf(extfile, "%llu %g %g %g %g %g\n", step, -2., yspan, zspan, shell_rgx2/num_shell_monos, mono_list[indentation2].get_prev_pos(0) - mono_list[indentation1].get_prev_pos(0));
else
fprintf(extfile, "%llu %g %g %g %g %g\n", step, -2., yspan, zspan, shell_rgx2/num_shell_monos, shell_extension);
}
else
{
fprintf(extfile, "%llu %g %g %g %g %g\n", step, 0., 0., 0., 0., 0.);
}

fflush(extfile);
}//if step>=0

//to make a movie, get rid of box resizing
//if(false)
if(BOX_RESIZE)
if(step % (NUMSKIP*100) == 0)
{
 resize_box(shell_extension, yspan, zspan);

 if(!THERMAL)
 {
  if((step > 15500000) && (fabs(shell_extension - saved_extension) < (0.05*SHELL_RADIUS)*1.e-6)) 
  {
	fprintf(stderr, "converged, exiting\n");

        char mvcommand[128];
        sprintf(mvcommand, "mv /tmp/[rs]*%6.6i output/", TRIALNUMBER);
        dummy= system(mvcommand);

	exit(1);
  }
  else
  {
   saved_extension = shell_extension;
  }

 }
}

#if PRINT_CM_ONLY
if(NUMBER_IN_POLYMER > 0)
{
fprintf(outfile, "%li %g %g %g\n", step, center_mass_pos[0], center_mass_pos[1], center_mass_pos[2]);

     fflush(outfile);
     fclose(outfile);
}
#endif


}


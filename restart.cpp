void write_restart(unsigned long long step)
{
int ii,kk;
FILE *restartfile;
char restartname[128];
int dummy;
if(RESTART)
 sprintf(restartname, "restart/restart_r_%6.6i", TRIALNUMBER);
else
 sprintf(restartname, "restart/restart_new_%6.6i", TRIALNUMBER);
 
restartfile = fopen(restartname, "w");
int loadmono;


fprintf(restartfile, "%g\n", pipette_dist_from_center);

for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
{
   fprintf(restartfile, "%i ", ii);
   for(kk = 0; kk < DIMENSION; kk++)
	fprintf(restartfile, "%g ", mono_list[ii].get_pos(kk));
   
   loadmono = (int) mono_list[ii].get_load_mono();
   fprintf(restartfile, "%i ", loadmono);
   if(loadmono == 1)
   {
      fprintf(restartfile, "%g ", mono_list[ii].get_sign_load());
   }
   else
      fprintf(restartfile, "%i ", 0);

   fprintf(restartfile, "\n");
}

for(ii = 0; ii < pairs.size(); ii++)
{
   fprintf(restartfile, "%i %li %li %g %i %i %g %g %g %g\n", ii, (*(pairs[ii].first)).get_id(), (*(pairs[ii].second)).get_id(), pairs[ii].get_bond_length(), (int) pairs[ii].get_branch_base(), (int) pairs[ii].get_parb_pair(), pairs[ii].get_bond_strength(), pairs[ii].get_bend_strength(), pairs[ii].get_orient_strength(), pairs[ii].get_nn_orient_strength());
}

fprintf(restartfile, "%llu\n", step);

fflush(restartfile);
fclose(restartfile);

}



//I haven't been recording load monos! I'll need to assume that load monos are still the closest to the poles... this should be okay, at least at high enough force...
void read_restart()
{
int ii, kk;
FILE *restartfile;
char restartname[128];
int dummy;
sprintf(restartname, "restart/restart_new_%6.6i", TRIALNUMBER);
restartfile = fopen(restartname, "r");

double temp_pos;
int temp_id;
int temp_indicator;
double temp_indicator_s;
int temp_id1, temp_id2;
double temp_length;
int temp_scrap_int;
double temp_strength;
double temp_scrap_double;

initialize_system_params();
initialize_files();
initialize_neighbor_list();

//fprintf(stderr, "???\n");
dummy= fscanf(restartfile, "%lg", &temp_pos);
//fprintf(stderr, "temppos %g\n", temp_pos);
pipette_dist_from_center = temp_pos;

double dcoeff;
for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
{
   dummy= fscanf(restartfile, "%i", &temp_id);
   //fprintf(stderr, "temp_id %i xyz", temp_id); 
   if(ii < NUMBER_IN_POLYMER)
	dcoeff= POLYMER_TDIFF_COEFF_FACTOR*TDIFF_COEFF;
   else
	dcoeff= TDIFF_COEFF;

   dynamic_monomer tempdmono(temp_id, dcoeff);
   mono_list.push_back(tempdmono);

   for(kk = 0; kk < DIMENSION; kk++)
   {
	dummy= fscanf(restartfile, "%lg", &temp_pos);
	//fprintf(stderr, " %g", temp_pos);
	mono_list[ii].set_pos(kk, temp_pos);
   }
   mono_list[ii].set_prev_pos();
   dummy= fscanf(restartfile, "%i %lg", &temp_indicator, &temp_indicator_s);
   //fprintf(stderr, "indicators %i %g\n", temp_indicator, temp_indicator_s);
   if(temp_indicator == 1)
	mono_list[ii].choose_load_mono(temp_indicator_s, pipette_dist_from_center);

  if(ii < NUMBER_IN_POLYMER)
  {
    mono_list[ii].set_exc_vol_spring_const(POLYMER_EXC_VOL_SPRING);
    mono_list[ii].set_radius(cl_polymer_mono_rad);
  }
}
//fprintf(stderr, "exited somehow\n");

while(fscanf(restartfile, "%i %i %i %lg %i %i %lg %lg %lg %lg", &temp_id, &temp_id1, &temp_id2, &temp_length, &temp_scrap_int, &temp_scrap_int, &temp_strength, &temp_scrap_double, &temp_scrap_double, &temp_scrap_double) != EOF) 
{
//fprintf(stderr, "%i %i %i %g %i %i %g %g %g %g\n", temp_id, temp_id1, temp_id2, temp_length, temp_scrap_int, temp_scrap_int, temp_strength, temp_scrap_double, temp_scrap_double, temp_scrap_double);
/*
if((temp_id1 > NUMBER_IN_POLYMER) && (temp_id2 > NUMBER_IN_POLYMER))
  temp_strength = choose_shell_spring_const(sqrt(mono_list[temp_id1].calculate_distance_sq(&mono_list[temp_id2], PREVIOUS)));*/

    monomer_pair temp_obj(&mono_list[temp_id1], &mono_list[temp_id2], NOT_BRANCH, false, temp_strength, 0., 0., 0.); 
//temp_strength is actually not really necessary because strength could be set by command line.. but it makes my life easier here.
    pairs.push_back(temp_obj);

}

//fprintf(stderr, "hi\n");
fflush(restartfile);
fclose(restartfile);

//write_state_to_pdb(0);
//exit(1);

update_system(0);

}



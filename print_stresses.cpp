void print_stresses(unsigned long long step)
{
	int ii,kk,jj,nn,mm, ll;
	int start_value;
	int mono1, mono2;
	double extension;
	double position[DIMENSION];
	double tot_tensile_strain;
        double tot_compressive_strain;
        vector<double> tot_compressive_energy(num_shell_monos, 0.);
        vector<double> tot_tensile_energy(num_shell_monos, 0.);
	vector<double> tot_elastic_energy(num_shell_monos, 0.);
        char tensile_strain_file_name[128];
	FILE *tensile_strain_file;

char cos_angle_file_name[128];
FILE *cos_angle_file;

double cos_angle;

int icounter, jcounter, kcounter;
int ix, jy, kz;
int ix2, jy2, kz2;
bool same_cell = false;
bool avoid = false;
double distance_between;

        char compressive_strain_file_name[128];
	FILE *compressive_strain_file;
        char tot_energy_file_name[128];
	FILE *tot_energy_file;

sprintf(tensile_strain_file_name, "output/stretch%6.6i_%llu", TRIALNUMBER, step);
tensile_strain_file = fopen(tensile_strain_file_name, "w");


sprintf(cos_angle_file_name, "output/cos_angles%6.6i_%llu",TRIALNUMBER, step);
cos_angle_file = fopen(cos_angle_file_name, "w");



        if(NUMBER_IN_POLYMER == 0)
          start_value = 0;
        else
          start_value = NUMBER_IN_POLYMER - 1;

	for(ii = start_value; ii < pairs.size(); ii++)
	{
          mono1 = (*(pairs[ii].first)).get_id();
          mono2 = (*(pairs[ii].second)).get_id();
          extension = pairs[ii].get_separation() - pairs[ii].get_bond_length();


          if(extension > 0.)
	  {
		tot_tensile_strain = extension;
		tot_elastic_energy[mono1-NUMBER_IN_POLYMER] += 0.25*SHELL_BOND_SPRING*extension*extension;//0.25 instead of 0.5 because I divide the energy contribution among 2 monomers
		tot_elastic_energy[mono2-NUMBER_IN_POLYMER] += 0.25*SHELL_BOND_SPRING*extension*extension;
		tot_tensile_energy[mono1-NUMBER_IN_POLYMER] += 0.25*SHELL_BOND_SPRING*extension*extension;
                tot_tensile_energy[mono2-NUMBER_IN_POLYMER] += 0.25*SHELL_BOND_SPRING*extension*extension;
	  }
	  else
	  {
		tot_tensile_strain = 0.;
	  }

	  for(kk = 0; kk < DIMENSION; kk++)
	    position[kk] = 0.5*(mono_list[mono1].get_prev_pos(kk)+mono_list[mono2].get_prev_pos(kk));

	  cos_angle = fabs(mono_list[mono1].get_prev_pos(0) - mono_list[mono2].get_prev_pos(0)) / extension;

	  fprintf(tensile_strain_file, "%g %g %g %g\n", position[0], position[1], position[2], tot_tensile_strain);
	  fprintf(cos_angle_file, "%g %g %g %g\n", position[0], position[1], position[2], cos_angle);
	}//loop over pairs for tension calc
fflush(tensile_strain_file);
fclose(tensile_strain_file);
fflush(cos_angle_file);
fclose(cos_angle_file);

sprintf(compressive_strain_file_name, "output/compression%6.6i_%llu", TRIALNUMBER, step);
compressive_strain_file = fopen(compressive_strain_file_name, "w");	

	
for(ii = 0; ii < N[0]; ii++)
{
 for(jj = 0; jj < N[1]; jj++)
 {
   for(kk = 0; kk < N[2]; kk++)
   {
    nn = get_firstmonoincell(ii,jj,kk);
    
    while((nn != -1) && (nn < NUMBER_IN_POLYMER))
	nn = get_monolinklist(nn); 

    while(nn != EMPTY)
    {
      icounter = 0;
      for(ix2 = (ii - 1 + N[0])%N[0]; icounter < 3; ix2 = (ix2+1)%N[0], icounter++)
      {
        ix = ix2;
        jcounter = 0;
        for(jy2 = (jj - 1 + N[1])%N[1]; jcounter < 3; jy2 = (jy2+1)%N[1], jcounter++)
        {
          jy = jy2;
          kcounter = 0;
          for(kz2 = kk; kcounter < 2; kz2 = (kz2+1)%N[2], kcounter++)
          {
            same_cell = false;
            kz = kz2;
            if(kz == kk)
            {
              if(jy2 == (jj-1+N[1])%N[1] || (jy == jj && ix2 == (ii-1+N[0])%N[0]))
              {
                continue;
              }//jy2 etc

              if(jy == jj)
               if(ix == ii)
                same_cell = true;
             }//kz==kk

            mm = get_firstmonoincell(ix,jy,kz);

            while((mm != -1) && (mm < NUMBER_IN_POLYMER))
              mm = get_monolinklist(mm);

            while(mm != LAST)
            {
              avoid = false;
              if(same_cell)
              if(mm <= nn)
               avoid = true;

              if(!avoid)
              {
		distance_between = (mono_list[nn]).calculate_distance_sq(&mono_list[mm], PREVIOUS);
                if(distance_between < MONO_DIAM2)
		{
		  tot_compressive_strain = sqrt(distance_between);
                  tot_elastic_energy[nn-NUMBER_IN_POLYMER] += 0.25*SPRINGCONST*tot_compressive_strain*tot_compressive_strain;//0.25 instead of 0.5 because I divide the energy contribution among 2 monomers
                  tot_elastic_energy[mm-NUMBER_IN_POLYMER] += 0.25*SPRINGCONST*tot_compressive_strain*tot_compressive_strain;
                  tot_compressive_energy[nn-NUMBER_IN_POLYMER] += 0.25*SPRINGCONST*tot_compressive_strain*tot_compressive_strain;
                  tot_compressive_energy[mm-NUMBER_IN_POLYMER] += 0.25*SPRINGCONST*tot_compressive_strain*tot_compressive_strain;

		  for(ll = 0; ll < DIMENSION; ll++)
		   position[ll] = 0.5*((mono_list[nn]).get_prev_pos(ll)+(mono_list[mm]).get_prev_pos(ll));
		  
		  fprintf(compressive_strain_file, "%g %g %g %g\n", position[0], position[1], position[2], tot_compressive_strain);
		}
		else
		{
		  tot_compressive_strain = 0.;
		}
              }//!avoid
	      do
	      {
               mm = get_monolinklist(mm);
	      }while((mm != -1) && (mm < NUMBER_IN_POLYMER));
            }//while mm!=LAST
          }//for(kz2=
        }//for(jy2=
      }//for(ix2=
      do
      {
        nn = get_monolinklist(nn);
      }while((nn != -1) && (nn < NUMBER_IN_POLYMER));
    }//while nn!= EMPTY
   }//for(kk
 }//for(jj
}//for(ii


fflush(compressive_strain_file);
fclose(compressive_strain_file);

sprintf(tot_energy_file_name, "output/energy%6.6i_%llu", TRIALNUMBER, step);
tot_energy_file = fopen(tot_energy_file_name, "w");


//need to print out mono specific total elastic energies.
for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
fprintf(tot_energy_file, "%g %g %g %g %g %g\n", mono_list[ii].get_prev_pos(0), mono_list[ii].get_prev_pos(1), mono_list[ii].get_prev_pos(2), tot_tensile_energy[ii-NUMBER_IN_POLYMER], tot_compressive_energy[ii-NUMBER_IN_POLYMER], tot_elastic_energy[ii-NUMBER_IN_POLYMER]);
}
fflush(tot_energy_file);
fclose(tot_energy_file);



}//end 





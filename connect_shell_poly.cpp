void connect_shell_poly()
{
int ii, jj, kk;
double temp;
double linkdist2 = 95.0*SHELL_RADIUS;
int linkmono;
double dist_from_center = 95.0*SHELL_RADIUS;

for(ii=0; ii < NUMBER_OF_SHELL_POLY_CONNECTIONS; ii++)
{
	dist_from_center = 95.0*SHELL_RADIUS;
	while(dist_from_center > (SHELL_RADIUS - 2.0)*(SHELL_RADIUS - 2.0))
	{	
          jj= (int)(unif_rand()*NUMBER_IN_POLYMER);
	  for(kk = 0; kk < DIMENSION; kk++)
	  {
	   dist_from_center = (mono_list[jj].get_prev_pos(kk)-L[kk]*0.5)*(mono_list[jj].get_prev_pos(kk)-L[kk]*0.5);
	  }
	}

        linkdist2 = 95.0*SHELL_RADIUS;

	if(!mono_list[jj].get_shell_poly_link())
	{
	  for(kk = NUMBER_IN_POLYMER; kk < NUMBER_OF_MONOMERS; kk++)
	  {
		if(!mono_list[kk].get_shell_poly_link())
		{
			temp = mono_list[jj].calculate_distance_sq(&mono_list[kk], PREVIOUS);
                        if(temp < linkdist2)
                        {
			  linkdist2 = temp;
                          linkmono = kk;
			}
		}//if kk !linked
	  }//for(kk

fprintf(stderr, "shell poly connect %i %i\n", jj, linkmono);
                monomer_pair temp_crosslinkpair(&mono_list[jj], &mono_list[linkmono], NOT_BRANCH, true, BOND_SPRING, 0., 0., 0.);
//I noticed this was bond length on 7/19/16... not changed yet though.... 
//                temp_crosslinkpair.set_bond_length(MONO_DIAM);
//prob should use this:
//changed on 8/16/16
		temp_crosslinkpair.set_bond_length(MONO_RAD+cl_polymer_mono_rad);
                crosslinkpairs.push_back(temp_crosslinkpair);
                mono_list[jj].set_shell_poly_link(true);
                mono_list[linkmono].set_shell_poly_link(true);

	}//if jj !crossli
	else
	{
		ii = ii - 1;
	}

}//for(ii







}


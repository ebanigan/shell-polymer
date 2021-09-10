/*
The purpose of this subroutine is to call the functions that update positions and
polarizations and update neighbor lists & pair objects
*/


void update_system(long tstep)
{

   long nn,ii,jj;
   int kk;
   double inv_numshell;
   double force;

	if(LOAD_TEST)
        if(spring_mono_id != -1)
	{ 
	  for(kk = 0; kk < DIMENSION; kk++)
	  {
	    force = -1.*LOAD_TEST_SPRING_CONST*(mono_list[spring_mono_id].get_prev_pos(kk) - spring_att_point[kk]);
	    force *= mono_list[spring_mono_id].get_tdiffusion_coeff()*invKT*dt;
	    mono_list[spring_mono_id].move_1d(kk, force);
	  }
	}

	 inv_numshell = 1./((double)num_shell_monos);
         for(kk = 0; kk < DIMENSION; kk++)
		shell_cm[kk] = 0.;

      for(nn = 0; nn < NUMBER_OF_MONOMERS; nn++)
      {
        (mono_list[nn]).check_translational_pbc();

        if(mono_list[nn].get_update_bool())
        {
 	 mono_list[nn].set_prev_pos();
//not using polarization in this code...
//	 mono_list[nn].set_prev_polarization();
        }
        else
	{
	   for(kk = 0; kk < DIMENSION; kk++)
		mono_list[nn].set_pos(kk, mono_list[nn].get_prev_pos(kk));
	}

	if(nn >= NUMBER_IN_POLYMER)
	 for(kk = 0; kk < DIMENSION; kk++)
	  shell_cm[kk] += mono_list[nn].get_prev_pos(kk);
      }
      for(kk = 0; kk < DIMENSION; kk++)
	shell_cm[kk] *= inv_numshell;


update_mono_list();



for(nn = 0; nn < pairs.size(); nn++)
  (pairs[nn]).update_relative_positions();
for(nn = 0; nn < crosslinkpairs.size(); nn++)
  crosslinkpairs[nn].update_relative_positions();



if(LENGTH_CONTROLLED_LOAD) // advance or relax the position of the pipette, test to see if condition to end simulation is met
{
  if(tstep > 0)
  {
    if(release_pipette(tstep))
    {
      pipette_dist_from_center -= dt*PIPETTE_VELOCITY;
      if(end_condition(tstep))
      {
          print_ext(tstep);
          exit(1);
      }
    }
    else
      pipette_dist_from_center += dt*PIPETTE_VELOCITY;
  }
  
}

if(COMPRESS)
{
BOTTOM_WALLPOS += BOTTOM_WALLFORCE + F_LOAD*dt*WALLDIFF;
TOP_WALLPOS += TOP_WALLFORCE - F_LOAD*dt*WALLDIFF;

TOP_WALLFORCE = 0.;
BOTTOM_WALLFORCE = 0.;
}

}




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

/*
 	if(num_shell_monos > 0)
	{
*/
	 inv_numshell = 1./((double)num_shell_monos);
         for(kk = 0; kk < DIMENSION; kk++)
		shell_cm[kk] = 0.;
/*	}
	else
	{
	inv_numshell = 0.;
         for(kk = 0; kk < DIMENSION; kk++)
                shell_cm[kk] = L[kk]*0.5;
	}
*/

      for(nn = 0; nn < NUMBER_OF_MONOMERS; nn++)
      {
        (mono_list[nn]).check_translational_pbc();

        if(mono_list[nn].get_update_bool())
        {
 	 mono_list[nn].set_prev_pos();
//don't use polarization in this code...
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


if(SOLID_INTERIOR)
{
	central_mono.set_prev_pos();
}

update_mono_list();



for(nn = 0; nn < pairs.size(); nn++)
  (pairs[nn]).update_relative_positions();
for(nn = 0; nn < crosslinkpairs.size(); nn++)
  crosslinkpairs[nn].update_relative_positions();

if(LENGTH_CONTROLLED_LOAD)
if(tstep > 0)
{
  pipette_dist_from_center += dt*PIPETTE_VELOCITY;
}

if(COMPRESS)
{
BOTTOM_WALLPOS += BOTTOM_WALLFORCE + F_LOAD*dt*WALLDIFF;
TOP_WALLPOS += TOP_WALLFORCE - F_LOAD*dt*WALLDIFF;

TOP_WALLFORCE = 0.;
BOTTOM_WALLFORCE = 0.;
}

//fprintf(stderr, "cm %g %g %g\n", shell_cm[0], shell_cm[1], shell_cm[2]);
}




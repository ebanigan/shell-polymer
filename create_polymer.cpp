void create_polymer()
{
     long ii,jj,ll;
     int kk;
     double initial_pos;
     double sum_sq;
     double random_displacement[DIMENSION];//used to get random walk for creation of polymer
     double prev_initial_pos[DIMENSION];
     double total[DIMENSION];
     double temp_pos;		

     for(kk = 0; kk < DIMENSION; kk++)
		prev_initial_pos[kk] = 0.5*L[kk];

	
//first monomer in polymer
    for(kk = 0; kk < DIMENSION; kk++)
    {
      (mono_list[0]).set_pos(kk, prev_initial_pos[kk]);
      (mono_list[0]).check_translational_pbc();//seems unnecessary?
    }
	(mono_list[0]).set_prev_pos();

if(NUMBER_IN_POLYMER > 1)
{
//mono_list 1 position is not yet set. not very elegant but practical solution is to just temporarily set mono1's position 
	for(kk = 0; kk < DIMENSION-1; kk++)
	  mono_list[1].set_pos(kk, mono_list[0].get_prev_pos(kk));
	kk=DIMENSION-1;//just to be sure..
	mono_list[1].set_pos(kk, mono_list[0].get_prev_pos(kk)+1.0);
	mono_list[1].set_prev_pos();

	(mono_list[0]).polymerize(&mono_list[1]);//so mono 0 should be the minus end of the polymer
	monomer_pair temp_pair(&mono_list[0], &mono_list[1], true, BOND_SPRING, POLYMER_STIFFNESS, POLYMER_POLARIZABILITY, POLYMER_NN_POLARIZATION_INT);
        (mono_list[1]).set_pair_num2(pairs.size());
	pairs.push_back(temp_pair);
}	 

//Main part of polymer
for(ii = 1; ii < NUMBER_IN_POLYMER; ii++)//it's okay to count up to NUMBER_IN_POLYMER since I have an if statement buried in the for loop
{
double next_pos[DIMENSION];
next_pos[0] = 9.;
next_pos[1] = 9.;
next_pos[2] = 9.*SHELL_RADIUS;

double xfactor, yfactor, zfactor;
if(!ELLIPSOID)
{
 xfactor = 1.;
 yfactor = 1.;
 zfactor = 1.;
}
else
{
 xfactor = (SHELL_RADIUS/MAJOR_AXIS)*(SHELL_RADIUS/MAJOR_AXIS);
 yfactor = (SHELL_RADIUS/MINOR_AXIS_1)*(SHELL_RADIUS/MINOR_AXIS_1);
 zfactor = (SHELL_RADIUS/MINOR_AXIS_2)*(SHELL_RADIUS/MINOR_AXIS_2);
}

while((next_pos[0]-0.5*LX)*(next_pos[0]-0.5*LX)*xfactor + (next_pos[1]-0.5*LY)*(next_pos[1]-0.5*LY)*yfactor + (next_pos[2]-0.5*LZ)*(next_pos[2]-0.5*LZ)*zfactor > SHELL_RADIUS*SHELL_RADIUS)
{
	   sum_sq = 0.;
	   for(kk = 0; kk < DIMENSION; kk++)
	   {
              random_displacement[kk] = 1.4*(unif_rand() - 0.5)*(2.0*cl_polymer_mono_rad);
              sum_sq += random_displacement[kk]*random_displacement[kk];
	   }
	   for(kk = 0; kk < DIMENSION; kk++)
           {
		random_displacement[kk] *= 2.00*cl_polymer_mono_rad/sqrt(3.*sum_sq);//why normalizing displacements? because the rest length between bonds is set by the initial separation
	        next_pos[kk] = prev_initial_pos[kk] + random_displacement[kk];
	   }
}//while
       for(kk = 0; kk < DIMENSION; kk++)
       {  
          //previous' position + distance of just barely touching monomers          
          initial_pos = next_pos[kk]; //prev_initial_pos[kk] + random_displacement[kk];
          (mono_list[ii]).set_pos(kk, initial_pos);
          (mono_list[ii]).check_translational_pbc();

          prev_initial_pos[kk] = initial_pos;
       }//for kk
	   (mono_list[ii]).set_prev_pos();
       //if monomer is not the last in the filament, polymerize the next one
       if(ii!=NUMBER_IN_POLYMER-1)
       {
        for(kk = 0; kk < DIMENSION-1; kk++)
            mono_list[ii+1].set_pos(kk, mono_list[ii].get_prev_pos(kk));
        mono_list[ii+1].set_pos(kk, mono_list[ii].get_prev_pos(kk)+1.0);
        mono_list[ii+1].set_prev_pos();

         (mono_list[ii]).polymerize(&mono_list[ii+1]);
         monomer_pair tempobj(&mono_list[ii], &mono_list[ii+1], true, BOND_SPRING, POLYMER_STIFFNESS, POLYMER_POLARIZABILITY, POLYMER_NN_POLARIZATION_INT);
	 (mono_list[ii+1]).set_pair_num2(pairs.size());
         pairs.push_back(tempobj);
       } // ii != NUMBER_IN_POLYMER-1
}//for(ii..

int num_cuts = (int)(POLYMER_CUT_PERCENTAGE*NUMBER_IN_POLYMER);
for(ii = 0; ii < num_cuts; ii++)
{
   jj = (int)(unif_rand()*pairs.size());
   fprintf(stderr, "cutting pair %li containing %li and %li\n", jj, (*(pairs[jj].first)).get_id(), (*(pairs[jj].second)).get_id());
   (*(pairs[jj].second)).set_pair_num2(-1);
   for(ll = jj + 1; ll < pairs.size(); ll++)
    (*(pairs[ll].second)).decrement_pair_num2();
   pairs.erase(pairs.begin()+jj);
}


if(NUMBER_IN_POLYMER > 0)
{
for(kk = 0; kk < DIMENSION; kk++)
{
  total[kk] = 0.;
}

for(ii = 0; ii < NUMBER_IN_POLYMER; ii++)
 for(kk = 0; kk < DIMENSION; kk++)
 {
  temp_pos = (mono_list[ii]).get_prev_pos(kk);
  total[kk] += temp_pos;
 }

for(kk = 0; kk < DIMENSION; kk++)
 total[kk] *= 1./(NUMBER_IN_POLYMER);

for(ii = 0; ii < NUMBER_IN_POLYMER; ii++)
{
  (mono_list[ii]).set_pos(0, (mono_list[ii]).get_pos(0) + LX*0.5 - total[0]);
  (mono_list[ii]).set_pos(1, (mono_list[ii]).get_pos(1) + LY*0.5 - total[1]);
  (mono_list[ii]).set_pos(2, (mono_list[ii]).get_pos(2) + LZ*0.5 - total[2]);
  (mono_list[ii]).set_prev_pos();
}
} 

double temp;
double linkdist2 = 95.0*SHELL_RADIUS;
int linkmono;
bool crosslink_created;

for(ii = 0; ii < NUMBER_OF_CROSSLINKS; ii++)
{
	jj= (int)(unif_rand()*NUMBER_IN_POLYMER);
	linkdist2 = 95.0*SHELL_RADIUS;

	if(!mono_list[jj].get_crosslinked())
	{
		crosslink_created = false;
		while(!crosslink_created)
		{
		 kk= (int)(unif_rand()*NUMBER_IN_POLYMER);
		 if(!mono_list[kk].get_crosslinked())
		 {
		  if((kk > jj+4) || (kk< jj-4))
		  {
			linkmono = kk;

		    fprintf(stderr, "creating crosslinkpair %li %i with bonds stiffness %g\n", jj, linkmono, cl_crosslink_spring_factor*BOND_SPRING);
		    monomer_pair temp_crosslinkpair(&mono_list[jj], &mono_list[linkmono], true, cl_crosslink_spring_factor*BOND_SPRING, 0., 0., 0.);
		    temp_crosslinkpair.set_bond_length(2.0*cl_polymer_mono_rad);
		    crosslinkpairs.push_back(temp_crosslinkpair);
		    mono_list[jj].set_crosslinked(true);
		    mono_list[linkmono].set_crosslinked(true);
		    crosslink_created = true;
		  }
		  else
			crosslink_created = false;//unnecessary, but do it anyway
		 }
		 else
		  crosslink_created = false;//unnecessary, but do it anyway
		}//while 

	}//if(!mono_list[jj].crosslinked()
	else
	{
		ii = ii - 1;
	}
}//for(ii


/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////

        double unit_vector[DIMENSION];
        double dist = 0.;
        double disp_vector[DIMENSION];
        double frac_str;

if(NUMBER_IN_POLYMER > 0)
for(ii = 0; ii < NUM_INIT_POLY_STEPS; ii++)
{
if(ii % (NUMSKIP*10) == 0)
  fprintf(stderr, "polymer initial dynamics step = %li\n", ii);

if(NUMBER_OF_CROSSLINKS > 0)
{
	frac_str = ((double)(ii+1.0))/((double)NUM_INIT_POLY_STEPS);
	relax_crosslinks(frac_str);//prevent huge jumps due to random crosslinking.
}

monomer_dynamics(true, 0);
polymer_interactions();

//idea is to relax polymer, while keeping it inside the shell.
for(jj = 0; jj < NUMBER_IN_POLYMER; jj++)
{
if((mono_list[jj].get_pos(0)-0.5*LX)*(mono_list[jj].get_pos(0)-0.5*LX) + (mono_list[jj].get_pos(1)-0.5*LY)*(mono_list[jj].get_pos(1)-0.5*LY) + (mono_list[jj].get_pos(2)-0.5*LZ)*(mono_list[jj].get_pos(2)-0.5*LZ) > (SHELL_RADIUS-0.5)*(SHELL_RADIUS-0.5))
{
	 dist = 0.;

	for(kk = 0; kk < DIMENSION; kk++)
	{
		unit_vector[kk] = mono_list[jj].get_pos(kk)-0.5*L[kk];
		dist += unit_vector[kk]*unit_vector[kk];
	}
	dist = sqrt(dist);
	for(kk = 0; kk < DIMENSION; kk++)
	{
		unit_vector[kk] *= 1./dist;
		disp_vector[kk] = unit_vector[kk] * (SHELL_RADIUS-0.5-dist);
	}
        mono_list[jj].move(disp_vector);
}
}//loop over jj<NUMBER_IN_POLYMER

update_system(ii);
}//loop over ii<numinit

}//end of subroutine



void relax_crosslinks(double fractional_str)
{
int ii;
for(ii = 0; ii < NUMBER_OF_CROSSLINKS; ii++)
     crosslinkpairs[ii].set_bond_strength(fractional_str*cl_crosslink_spring_factor*BOND_SPRING);
}


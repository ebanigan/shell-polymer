void create_cylinder()
{
int ii, jj, kk;

double rand_angle;
double rand_pos;

double cyl_area = TWOPI*SHELL_RADIUS*CYLINDER_LENGTH;
double bond_length = sqrt(cyl_area / num_shell_monos);

double increment_angle = 2.*asin(bond_length / SHELL_RADIUS);
num_boundary =(int) (TWOPI / increment_angle);
increment_angle = TWOPI / (double) num_boundary;

cl_num_load_monos = num_boundary;

double norm, sumsq;
double init_pos[DIMENSION];
init_pos[0] = 0.;

vector<bool> tempvec(num_shell_monos,false);
double bonddist, prev_bonddist;
double shellspringconst;
double low_bound;

for(ii = NUMBER_IN_POLYMER; ii < NUMBER_IN_POLYMER + num_boundary; ii++)
{
	jj = ii - NUMBER_IN_POLYMER;
	mono_list[ii].set_pos(0, 0.5*L[0] - 0.5*CYLINDER_LENGTH);
	mono_list[ii].set_pos(1, 0.5*L[1] + SHELL_RADIUS*cos(jj*increment_angle));
	mono_list[ii].set_pos(2, 0.5*L[2] + SHELL_RADIUS*sin(jj*increment_angle));
	mono_list[ii].set_prev_pos();
}

for(ii = NUMBER_IN_POLYMER + num_boundary; ii < NUMBER_IN_POLYMER + 2*num_boundary; ii++)
{
        jj = ii - NUMBER_IN_POLYMER - num_boundary;
        mono_list[ii].set_pos(0, 0.5*L[0] + 0.5*CYLINDER_LENGTH);
        mono_list[ii].set_pos(1, 0.5*L[1] + SHELL_RADIUS*cos(jj*increment_angle));
        mono_list[ii].set_pos(2, 0.5*L[2] + SHELL_RADIUS*sin(jj*increment_angle));
        mono_list[ii].set_prev_pos();
}

for(ii = NUMBER_IN_POLYMER + 2*num_boundary; ii < NUMBER_OF_MONOMERS; ii++)
{
	rand_angle = unif_rand()*TWOPI;
	rand_pos = 0.9999*unif_rand()*CYLINDER_LENGTH;
        mono_list[ii].set_pos(0, rand_pos - 0.5*CYLINDER_LENGTH + 0.5*L[0]);
	mono_list[ii].set_pos(1, SHELL_RADIUS*cos(rand_angle) + 0.5*L[1]);
	mono_list[ii].set_pos(2, SHELL_RADIUS*sin(rand_angle) + 0.5*L[2]);
        mono_list[ii].set_prev_pos();
}//for loop ii over  monos

update_mono_list();

fprintf(stderr, "beginning shell relaxation %lu\n", pairs.size());
for(ii = 0; ii < 2*NUM_RELAXATION_STEPS; ii++)
{
//write_state_to_pdb(ii);
if(ii % NUMSKIP == 0)
        fprintf(stderr, "relaxation %i\n", ii);
if(!SPRINGS_ONLY)
{
        monomer_dynamics(true,0);//let monomers relax a bit
        filament_interactions();//to keep the polymer in a string

	for(jj = NUMBER_IN_POLYMER; jj < NUMBER_IN_POLYMER+2*num_boundary; jj++)
	{
		for(kk = 0; kk < DIMENSION; kk++)
		  mono_list[jj].set_pos(kk, mono_list[jj].get_prev_pos(kk));
	}//keep boundary monos fixed
}

//keep other monos on cylinder surface
        for(jj = NUMBER_IN_POLYMER+2*num_boundary; jj < NUMBER_OF_MONOMERS; jj++)
        {
		init_pos[1] = mono_list[jj].get_pos(1) - 0.5*L[1];
                init_pos[2] = mono_list[jj].get_pos(2) - 0.5*L[2];
		sumsq = init_pos[1]*init_pos[1] + init_pos[2]*init_pos[2];
		norm = sqrt(sumsq);
		init_pos[1] *= SHELL_RADIUS / norm;
		init_pos[2] *= SHELL_RADIUS / norm;
		mono_list[jj].set_pos(1, init_pos[1] + 0.5*L[1]);
		mono_list[jj].set_pos(2, init_pos[2] + 0.5*L[2]);
		if(mono_list[jj].get_pos(0) > 0.5*(L[0] + CYLINDER_LENGTH))
		  mono_list[jj].set_pos(0, 0.5*(L[0] + CYLINDER_LENGTH) - MONO_DIAM);
		else if(mono_list[jj].get_pos(0) < 0.5*(L[0] - CYLINDER_LENGTH))
		  mono_list[jj].set_pos(0, 0.5*(L[0] - CYLINDER_LENGTH) + MONO_DIAM);
		mono_list[jj].set_prev_pos();
	}//jj loop

update_system(ii);
}//end for loop over relaxation steps

fprintf(stderr, "initial relaxation complete\n");


///*********************************************************************
//******************************************************************************************************
//************************************************************************************//
for(int aa = 0; aa < num_shell_monos; aa++)
{
        attachment_matrix.push_back(tempvec);
}
for(ii = 0; ii < num_shell_monos; ii++)
{
 num_attachments.push_back(0);
 for(jj = 0; jj < num_shell_monos; jj++)
        attachment_matrix[ii][jj] = false;
}

tot_attachments = 0;

for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
if(!(SPRINGS_ONLY || VARIABLE_BOND_LENGTH))
 bonddist = 0.9;
else
 bonddist = 0.75*sqrt(4.*PI*SHELL_RADIUS*SHELL_RADIUS / num_shell_monos);

prev_bonddist = 0.;

if(ii < NUMBER_IN_POLYMER + 2*num_boundary)
 low_bound = 3;
else
 low_bound = LOWER_BOUND_CONNECTIVITY;

while(num_attachments[ii-NUMBER_IN_POLYMER] < low_bound)//if fewer than lower bound connections, connect nearby monos.
{
for(jj = NUMBER_IN_POLYMER; jj < NUMBER_OF_MONOMERS; jj++)
{
if(mono_list[ii].calculate_distance_sq(&mono_list[jj], PREVIOUS) < bonddist*bonddist)
{
if(mono_list[ii].calculate_distance_sq(&mono_list[jj], PREVIOUS) > prev_bonddist*prev_bonddist)
{
if(!attachment_matrix[ii-NUMBER_IN_POLYMER][jj-NUMBER_IN_POLYMER])
{
if(num_attachments[ii-NUMBER_IN_POLYMER] < UPPER_BOUND_CONNECTIVITY)//don't let ii go above the upper bound... 
if(num_attachments[jj-NUMBER_IN_POLYMER] < UPPER_BOUND_CONNECTIVITY)//also check that jj is below upper bound...
{
//the reason for printing this way is easy mathematica visualization
      fprintf(stderr, "{%i,%i},", ii+1-NUMBER_IN_POLYMER, jj+1-NUMBER_IN_POLYMER);
        num_attachments[ii-NUMBER_IN_POLYMER] = num_attachments[ii-NUMBER_IN_POLYMER] + 1;
        num_attachments[jj-NUMBER_IN_POLYMER] = num_attachments[jj-NUMBER_IN_POLYMER] + 1;
        attachment_matrix[ii-NUMBER_IN_POLYMER][jj-NUMBER_IN_POLYMER] = true;
        attachment_matrix[jj-NUMBER_IN_POLYMER][ii-NUMBER_IN_POLYMER] = true;

        tot_attachments++;

        //pairing.
        shellspringconst = choose_shell_spring_const(sqrt(mono_list[ii].calculate_distance_sq(&mono_list[jj], PREVIOUS)));

        monomer_pair temp_obj(&mono_list[ii], &mono_list[jj], NOT_BRANCH, false, shellspringconst, 0., 0., 0.);
        //mono_list[jj].set_pair_num2(pairs.size()); // don't bother because many monos will be in multiple pairs

        pairs.push_back(temp_obj);


}//force num attached to be less than UPPER_BOUND_CONNECTIVITY
}//if not already attached
}//if dist is right
}

}//for(jj

//fprintf(stderr, "%i, %g\n", ii, bonddist);
//
/*if(bonddist > SHELL_RADIUS)
 *   fprintf(stderr, "bonddist larger than SHELL_RADIUS = %g\n", SHELL_RADIUS);
 *   else if(bonddist > 2.00)
 *     fprintf(stderr, "bonddist is now larger than 2\n");
 *     */


prev_bonddist = bonddist;
if(!(SPRINGS_ONLY || VARIABLE_BOND_LENGTH))
 bonddist = bonddist + 0.1;
else
 bonddist = bonddist + 0.075*sqrt(4.*PI*SHELL_RADIUS*SHELL_RADIUS / num_shell_monos);
if(bonddist > 10.)
{
  //fprintf(stderr, "bonddist for shell is larger than 10!\n");
  //// exit(1);
}

}//while
}//for(ii
fprintf(stderr, "\n");


/*Now check to see that all pieces of shell are attached*/
bool connected = false;
while(!connected)
{
fprintf(stderr, "checking connectedness\n");
connected = true;

for(ii = 0; ii < num_shell_monos; ii++)
 main_shell.push_back(false);

if(main_shell.size() > 0)
{
main_shell[0] = true;
visit(0);
}
else
  connected = true;

int new_bonds_to_add = 2;
int new_bonds_added = 0;

for(ii = 0; ii < num_shell_monos; ii++)
{
  if(!main_shell[ii])
  {
    fprintf(stderr, "disconn %i\n", ii);

    connected = false;
    //add a bond to the nearest non-maxed shell mono
    double min_dist2 = 100.*SHELL_RADIUS*SHELL_RADIUS;
    double min_mono = 0;

    if(num_attachments[ii] < UPPER_BOUND_CONNECTIVITY)
    for(jj = 0; jj < num_shell_monos; jj++)
    {
        if(!attachment_matrix[ii][jj])
        if(main_shell[jj])
        if(num_attachments[ii] < UPPER_BOUND_CONNECTIVITY)
        if(mono_list[ii+NUMBER_IN_POLYMER].calculate_distance_sq(&mono_list[jj+NUMBER_IN_POLYMER], PREVIOUS) < min_dist2)
        {
                min_mono = jj;
                min_dist2 = mono_list[ii+NUMBER_IN_POLYMER].calculate_distance_sq(&mono_list[jj+NUMBER_IN_POLYMER], PREVIOUS);
        }
    }//for(jj...


    add_bond(ii,min_mono);
    new_bonds_added++;

    if(new_bonds_added >= new_bonds_to_add)
    {
      main_shell.clear();
      break;
    }
  }//if(!main_shell[ii]) 
}//for(ii...

}//while shell is not connected.


fprintf(stderr, "tot attachments %i tot monos %i\n", tot_attachments, num_shell_monos);

attachment_matrix.clear();
main_shell.clear();

fprintf(stderr, "end create shell\n");


}//end create_cylinder


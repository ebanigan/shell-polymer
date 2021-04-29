void create_shell()
{
int ii,jj,kk,ii2;
double sum_sq = 99.;
double norm;
vector<double> init_pos(DIMENSION, 0.);
vector<bool> tempvec(num_shell_monos,false);
double bonddist, prev_bonddist;
double shellspringconst;


if(DIMENSION==3)
{

for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
if(ii < NUMBER_IN_POLYMER+NUM_BINDING_MONOS)
	mono_list[ii].choose_binding_mono();

sum_sq = 99.;
while(sum_sq > 1.)
{
  sum_sq = 0.;
  for(kk = 0; kk < DIMENSION; kk++)
  {
	init_pos[kk] = 2.*(unif_rand()-0.5);
	sum_sq = sum_sq + init_pos[kk]*init_pos[kk];
  }
}

norm = 2./sqrt(sum_sq)*SHELL_RADIUS;//starting things out at twice the desired shell radius
for(kk = 0; kk < DIMENSION; kk++)
{
	init_pos[kk] = init_pos[kk]*norm;
	init_pos[kk] = init_pos[kk]+L[kk]*0.5;
	//fprintf(stderr, "%i %i %g\n", ii, kk, init_pos[kk]);
	mono_list[ii].set_pos(kk, init_pos[kk]);
}
mono_list[ii].set_prev_pos();
}
}//shell set up for dimension=3
else if(DIMENSION==2)
{
}//dimension=2
else
 fprintf(stderr, "weird dimension. exit\n");

update_mono_list();

fprintf(stderr, "beginning shell relaxation %i\n", int(pairs.size()));

int relaxation_factor  = 1;
if(EXTENDED_RELAXATION)
  relaxation_factor = 10;

for(ii = 0; ii < relaxation_factor*NUM_RELAXATION_STEPS; ii++)
{
//write_state_to_pdb(ii);
if(ii % (4*NUMSKIP) == 0)
	fprintf(stderr, "relaxation %i\n", ii);

if((!SPRINGS_ONLY) || EXTENDED_RELAXATION)
{
	monomer_dynamics(true,0);//let monomers relax a bit
	filament_interactions();//to keep the polymer in a string
	update_system(ii);
}

	for(jj = NUMBER_IN_POLYMER; jj < NUMBER_OF_MONOMERS; jj++)
	{
		sum_sq = 0.;
		for(kk = 0; kk < DIMENSION; kk++)
		{
			init_pos[kk] = mono_list[jj].get_pos(kk) - 0.5*L[kk];
			sum_sq = sum_sq + (mono_list[jj].get_pos(kk) - 0.5*L[kk])*(mono_list[jj].get_pos(kk) - 0.5*L[kk]);
		}

		if(ii < relaxation_factor*NUM_RELAXATION_STEPS / 2)
		{
		   norm = (2.0-(double)ii/(relaxation_factor*NUM_RELAXATION_STEPS/2))*1./sqrt(sum_sq)*SHELL_RADIUS;
		}
		else
		{
		  norm = 1./sqrt(sum_sq)*SHELL_RADIUS;
		}

		for(kk = 0; kk < DIMENSION; kk++)
		{
        		init_pos[kk] = init_pos[kk]*norm;
        		init_pos[kk] = init_pos[kk]+L[kk]*0.5;
        		mono_list[jj].set_pos(kk, init_pos[kk]);
		}//put them back on the sphere

		mono_list[jj].set_prev_pos();
	}

for(jj = 0; jj < NUMBER_IN_POLYMER; jj++)//yes, this is in fact the correct cycle because I set shell prev pos above.
	mono_list[jj].set_prev_pos();

update_mono_list();
//write_state_to_pdb(ii);

if(ii % (25*NUMSKIP) == 0)
  run_resize_box();

}
fprintf(stderr, "initial relaxation complete\n");

choose_shell_attachments();

}//end create_shell()



void visit(int site_number)
{
 int jj;
 for(jj = 0; jj < num_shell_monos; jj++)
 {
   if(!main_shell[jj])
   if(attachment_matrix[site_number][jj])
   {
//     fprintf(stderr, "%i\n", jj);
     main_shell[jj] = true;
     visit(jj);
   }
 }
}


void add_bond(int ii, int jj)
{
int kk;
double shellspringconst;

     fprintf(stderr, "{%i,%i},", ii+1, jj+1);
 num_attachments[ii] = num_attachments[ii] + 1;
 num_attachments[ii] = num_attachments[jj] + 1;

        attachment_matrix[ii][jj] = true;
        attachment_matrix[jj][ii] = true;

        tot_attachments++;

        if(!LENGTH_DEPENDENT_SPRINGS)
                shellspringconst = SHELL_BOND_SPRING;
        else
        {
                shellspringconst = SHELL_BOND_SPRING / sqrt(mono_list[ii+NUMBER_IN_POLYMER].calculate_distance_sq(&mono_list[jj+NUMBER_IN_POLYMER], PREVIOUS));
                if(shellspringconst > LDEP_FACTOR*SHELL_BOND_SPRING)
	                  shellspringconst = LDEP_FACTOR*SHELL_BOND_SPRING;
                  else if(shellspringconst < SHELL_BOND_SPRING/LDEP_FACTOR)
                          shellspringconst = SHELL_BOND_SPRING/LDEP_FACTOR;
        }

        //pairing.
        monomer_pair temp_obj(&mono_list[ii+NUMBER_IN_POLYMER], &mono_list[jj+NUMBER_IN_POLYMER], false, shellspringconst, 0., 0., 0.);
        pairs.push_back(temp_obj);
}


double choose_shell_spring_const(double distance_between)
{
double shellspringconst;
        if(!LENGTH_DEPENDENT_SPRINGS)
                shellspringconst = SHELL_BOND_SPRING;
        else
        {
                shellspringconst = SHELL_BOND_SPRING / distance_between;
                if(shellspringconst > LDEP_FACTOR*SHELL_BOND_SPRING)
                        shellspringconst = LDEP_FACTOR*SHELL_BOND_SPRING;
                else if(shellspringconst < SHELL_BOND_SPRING/LDEP_FACTOR)
                        shellspringconst = SHELL_BOND_SPRING/LDEP_FACTOR;
        }
return shellspringconst;
}



void choose_shell_attachments()
{
int aa, ii, jj;
vector<bool> tempvec(num_shell_monos,false);
double bonddist, prev_bonddist;
double shellspringconst;

for(aa = 0; aa < num_shell_monos; aa++)
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
while(num_attachments[ii-NUMBER_IN_POLYMER] < LOWER_BOUND_CONNECTIVITY)//if fewer than lower bound connections, connect nearby monos.
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
        num_attachments[ii-NUMBER_IN_POLYMER] = num_attachments[ii-NUMBER_IN_POLYMER] + 1;
        num_attachments[jj-NUMBER_IN_POLYMER] = num_attachments[jj-NUMBER_IN_POLYMER] + 1;
        attachment_matrix[ii-NUMBER_IN_POLYMER][jj-NUMBER_IN_POLYMER] = true;
        attachment_matrix[jj-NUMBER_IN_POLYMER][ii-NUMBER_IN_POLYMER] = true;

        tot_attachments++;

        //pairing.
        shellspringconst = choose_shell_spring_const(sqrt(mono_list[ii].calculate_distance_sq(&mono_list[jj], PREVIOUS)));

	monomer_pair temp_obj(&mono_list[ii], &mono_list[jj], false, shellspringconst, 0., 0., 0.);
	pairs.push_back(temp_obj);
}//force num attached to be less than UPPER_BOUND_CONNECTIVITY
}//if not already attached
}//if dist is right
}//if dist is right, again
}

prev_bonddist = bonddist;
if(!(SPRINGS_ONLY || VARIABLE_BOND_LENGTH))
 bonddist = bonddist + 0.1;
else
 bonddist = bonddist + 0.075*sqrt(4.*PI*SHELL_RADIUS*SHELL_RADIUS / num_shell_monos);
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

int new_bonds_to_add = 3;
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

fprintf(stderr, "end choose shell attachments\n");

}//end choose_shell_attachments()



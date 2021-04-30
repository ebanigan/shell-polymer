void create_ellipsoid()
{
int ii, jj, kk;
double axes[DIMENSION];
double axes2[DIMENSION];
double axes4[DIMENSION];
double axes6[DIMENSION];
double positions[DIMENSION];
double sum_sq;
double inv_prob_max =1./(MAJOR_AXIS * MINOR_AXIS_1);
double prob;
bool position_accepted;
double init_pos[DIMENSION];
double init_pos2[DIMENSION];
double direction_vector[DIMENSION];
double sum_2nd, sum_4th, sum_6th;
double mult_factor;


axes[0] = MAJOR_AXIS;
axes[1] = MINOR_AXIS_1;
axes[2] = MINOR_AXIS_2;
for(kk = 0; kk < DIMENSION; kk++)
{
	axes2[kk] = axes[kk]*axes[kk];
	axes4[kk] = axes2[kk]*axes2[kk];
	axes6[kk] = axes4[kk]*axes2[kk];
}

for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
position_accepted = false;
while(!position_accepted)
{
sum_sq = 0.;
for(kk = 0; kk < DIMENSION; kk++)
{
 positions[kk]=100.;
 sum_sq += positions[kk]*positions[kk];
}
while(sum_sq > 1.)
{
sum_sq = 0.;
for(kk = 0; kk < DIMENSION; kk++)
{
  positions[kk]=2.*(unif_rand()-0.5);
  sum_sq += positions[kk]*positions[kk];
}
}//while

prob = inv_prob_max * sqrt( (MAJOR_AXIS*MINOR_AXIS_2*positions[1])*(MAJOR_AXIS*MINOR_AXIS_2*positions[1]) + (MAJOR_AXIS*MINOR_AXIS_1*positions[2])*(MAJOR_AXIS*MINOR_AXIS_1*positions[2]) + (MINOR_AXIS_1*MINOR_AXIS_2*positions[0])*(MINOR_AXIS_1*MINOR_AXIS_2*positions[0]) );


if(unif_rand() < prob)
{
for(kk = 0; kk < DIMENSION; kk++)
{
  positions[kk] *= axes[kk];
  mono_list[ii].set_pos(kk, 0.5*L[kk] + positions[kk]);
}
mono_list[ii].set_prev_pos();

position_accepted = true;
}
}//while(!position_accepted)

}//for(ii...


update_mono_list();

fprintf(stderr, "beginning ellipsoid shell relaxation %i\n", int(pairs.size()));

int relaxation_factor  = 1;
if(EXTENDED_RELAXATION)
  relaxation_factor = 10;


for(ii = 0; ii < relaxation_factor*NUM_RELAXATION_STEPS; ii++)
{
if(ii % NUMSKIP == 0)
        fprintf(stderr, "relaxation %i\n", ii);

if((!SPRINGS_ONLY) || EXTENDED_RELAXATION)
{
        monomer_dynamics(true,0);//let monomers relax a bit
        polymer_interactions();//to keep the polymer in a string
        update_system(ii);
}

//instead of having the shell shrink as relaxation proceeds, let have the excluded volume spring constant gradually increase
//actually, I don't even see where I do this. loks like excluded volume is constant.
        for(jj = NUMBER_IN_POLYMER; jj < NUMBER_OF_MONOMERS; jj++)
        {
//                sum_sq = 0.;
		sum_2nd = 0.;
		sum_4th = 0.;
		sum_6th = 0.;
                for(kk = 0; kk < DIMENSION; kk++)
                {
                        init_pos[kk] = mono_list[jj].get_pos(kk) - 0.5*L[kk];
			init_pos2[kk] = init_pos[kk]*init_pos[kk];
			direction_vector[kk] = init_pos[kk] / axes2[kk];
			sum_2nd += init_pos2[kk] / axes2[kk];
			sum_4th += init_pos2[kk] / axes4[kk];
			sum_6th += init_pos2[kk] / axes6[kk];
//                        sum_sq = sum_sq + (mono_list[jj].get_pos(kk) - 0.5*L[kk])*(mono_list[jj].get_pos(kk) - 0.5*L[kk]) / (axes[kk]*axes[kk]);
                }//maybe want to make sure denom > 0
		
		if(sum_6th > 1.e-12)
		{

//this code takes monomer from its position on a different ellipsoid (that it has fluctuated to) to a position on the ellipsoid I want
		 mult_factor = (sum_4th - sqrt(sum_4th*sum_4th - sum_6th*(-1.+sum_2nd))) / sum_6th; 
		 for(kk = 0; kk < DIMENSION; kk++)
		 {
		    init_pos[kk] -= mult_factor*direction_vector[kk];
		    init_pos[kk] += 0.5*L[kk];
		    mono_list[jj].set_pos(kk, init_pos[kk]);
		 }
		}
                 mono_list[jj].set_prev_pos();
	}//for loop over shell monos.

	for(jj = 0; jj < NUMBER_IN_POLYMER; jj++)
		mono_list[jj].set_prev_pos();

update_mono_list();

}//end of relaxation loop

fprintf(stderr, "initial relaxation complete\n");

choose_shell_attachments();
}//end create_ellipsoid



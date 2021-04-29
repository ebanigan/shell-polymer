void binding(int mono1, int mono2)
{
	double distance[DIMENSION];
	double force1[DIMENSION];
	double force2[DIMENSION];
	double distancesq = mono_list[mono1].calculate_distance_sq(&mono_list[mono2], true);
	int kk;
	double prefactor;// = TETHER_STR;

//was designing a function for non-specific binding
//abandoned
//now just a placeholder

/*
if(distancesq < TETHER_RANGE_SQ)
{
if(distancesq > 2.*cl_polymer_mono_rad)//MONO_DIAM)
{

	for(kk = 0; kk < DIMENSION; kk++)
	{
		distance[kk] = mono_list[mono1].calculate_1d_sep(&mono_list[mono2], kk, true);// pos mono 1 minus pos mono 2
		force1[kk] = -1.*(mono_list[mono1].get_tdiffusion_coeff())*dt*TETHER_STR*(-14.*distance[kk] * LJ_MIN14/pow(distancesq, 8) + 12.*distance[kk]*LJ_MIN12 / pow(distancesq, 7));//due to distance convention above, no neg. sign in front of force on this line, since force is force on mono1
		force2[kk] = -1.*force1[kk];
	}


	mono_list[mono1].move(force1);
	mono_list[mono2].move(force2);	

}
}
*/

}


void binding_null(int mono1, int mono2){}


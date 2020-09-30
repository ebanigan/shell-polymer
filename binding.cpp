void binding(int mono1, int mono2)
{
	double distance[DIMENSION];
	double force1[DIMENSION];
	double force2[DIMENSION];
	double distancesq = mono_list[mono1].calculate_distance_sq(&mono_list[mono2], true);
	int kk;
	double prefactor = TETHER_STR;

//fprintf(stderr, "in binding\n");

//if(EPSILON > 1e-5)
//{
if(distancesq < TETHER_RANGE_SQ)
{
if(distancesq > 2.*cl_polymer_mono_rad)//MONO_DIAM)
{
if(true)//if((mono1 != mono2 -1) && (mono1 != mono2+1))
{

//fprintf(stderr, "calculating binding, %i %i\n", mono1, mono2);

	for(kk = 0; kk < DIMENSION; kk++)
	{
		distance[kk] = mono_list[mono1].calculate_1d_sep(&mono_list[mono2], kk, true);// pos mono 1 minus pos mono 2
		force1[kk] = -1.*(mono_list[mono1].get_tdiffusion_coeff())*dt*TETHER_STR*(-14.*distance[kk] * LJ_MIN14/pow(distancesq, 8) + 12.*distance[kk]*LJ_MIN12 / pow(distancesq, 7));//due to distance convention above, no neg. sign in front of force on this line, since force is force on mono1
		force2[kk] = -1.*force1[kk];
	}


	mono_list[mono1].move(force1);
	mono_list[mono2].move(force2);	

/*
fprintf(stderr, "particle 1 (%g, %g, %g)\n", mono_list[mono1].get_prev_pos(0), mono_list[mono1].get_prev_pos(1), mono_list[mono1].get_prev_pos(2));
fprintf(stderr, "particle 2 (%g, %g, %g)\n", mono_list[mono2].get_prev_pos(0), mono_list[mono2].get_prev_pos(1), mono_list[mono2].get_prev_pos(2));
fprintf(stderr, "distancesq %g distance %g\n", distancesq, sqrt(distancesq));
fprintf(stderr, "force 1 (%g, %g, %g)\n", force1[0], force1[1], force1[2]);
fprintf(stderr, "force 2 (%g, %g, %g)\n", force2[0], force2[1], force2[2]);
int wait1;
scanf("%i", &wait1);
*/
/*exit(1);*/
}
}
}
//}

}


void binding_null(int mono1, int mono2){}


void print_shell_spring_data(unsigned long long step)
{
	//there are NUMBER_IN_POLYMER-1 springs connecting the monomers of the linear polymer
	//first element of pairs to pay attention to: pairs[NUMBER_IN_POLYMER]
	int ii;
	int start_value;
	int mono1, mono2;
	double extension;
	double cos_angle_with_force;
	double total_extension = 0.;
	double total_cos_bond_angle = 0.;//Note: really, abs cos bond angle, since 0 and 180 deg are equiv. for my purposes
	double total_parallel_extension = 0.;
	int num_bonds = 0;
	double avg_extension;
	double avg_cos_bond_angle;
	double avg_parallel_extension;

	if(NUMBER_IN_POLYMER == 0)
	  start_value = 0;
	else 
	  start_value = NUMBER_IN_POLYMER - 1;

//calculate cosine of spring angle with force
//calculate average extension of springs
//calculate average extension parallel to force axis
	for(ii = start_value; ii < pairs.size(); ii++)
	{

		mono1 = (*(pairs[ii].first)).get_id();
		mono2 = (*(pairs[ii].second)).get_id();

		if(!(mono_list[mono1].get_load_mono() && mono_list[mono2].get_load_mono()))
		{
		  cos_angle_with_force = fabs(pairs[ii].get_relative_position(0) * pairs[ii].get_inv_separation());
		  extension = pairs[ii].get_separation() - pairs[ii].get_bond_length();

		  total_extension = total_extension + extension;
		  total_cos_bond_angle = total_cos_bond_angle + cos_angle_with_force; 
		  total_parallel_extension = total_parallel_extension + extension*cos_angle_with_force;
		  num_bonds = num_bonds + 1;
		}
	}
	
	avg_extension = total_extension / ((double) num_bonds);
	avg_cos_bond_angle = total_cos_bond_angle / ((double) num_bonds);
	avg_parallel_extension = total_parallel_extension / ((double) num_bonds);

	fprintf(springfile, "%llu %g %g %g\n", step, avg_extension, avg_cos_bond_angle, avg_parallel_extension);
	fflush(springfile);


}




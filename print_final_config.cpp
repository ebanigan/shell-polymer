void print_final_config()
{
int ii, kk;
int mono1, mono2;
char configname[72];
FILE *configfile;
int num_poly_pairs;
double init_sep[DIMENSION];
double init_pos[DIMENSION];
double final_sep[DIMENSION];
double final_pos[DIMENSION];

if(NUMBER_IN_POLYMER < 2)
   num_poly_pairs = 0;
else
   num_poly_pairs = NUMBER_IN_POLYMER-1;

sprintf(configname, "output/config%6.6i", TRIALNUMBER);
configfile = fopen(configname, "w");

for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
	fprintf(configfile, "%i %g %g %g %g %g %g\n", ii, mono_list[ii].get_x0(0), mono_list[ii].get_x0(1), mono_list[ii].get_x0(2), mono_list[ii].get_pos(0), mono_list[ii].get_pos(1), mono_list[ii].get_pos(2));
}

for(ii = num_poly_pairs; ii < pairs.size(); ii++) 
{
	mono1 = (*(pairs[ii].first)).get_id();
        mono2 = (*(pairs[ii].second)).get_id();
	
	for(kk = 0; kk < DIMENSION; kk++)
	{
	  init_sep[kk] = fabs(mono_list[mono1].get_x0(kk) - mono_list[mono2].get_x0(kk));
	  init_pos[kk] = 0.5*(mono_list[mono1].get_x0(kk) + mono_list[mono2].get_x0(kk));
	  final_sep[kk] = fabs(mono_list[mono1].get_pos(kk) - mono_list[mono2].get_pos(kk));
	  final_pos[kk] = 0.5*(mono_list[mono1].get_pos(kk) + mono_list[mono2].get_pos(kk));
	} 

	fprintf(configfile, "%i %i %g %g %g %g %g %g %g %g %g %g %g %g\n", mono1, mono2, init_pos[0], init_pos[1], init_pos[2], final_pos[0], final_pos[1], final_pos[2], init_sep[0], init_sep[1], init_sep[2], final_sep[0], final_sep[1], final_sep[2]);
}
fflush(configfile);
fclose(configfile);

}


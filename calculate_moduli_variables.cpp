/*
-turn off size changes after initial transients -- do this from command line
-periodically reset x0 -- for now i won't do this
-run it, write code to calculate two moduli
-maybe a bool to turn on and off bend calc
*/

void calculate_moduli_variables(unsigned long long step)
{
int ii, kk;

/*
Calculate <R^2> and <R> at this moment (avg. over particles), where R is the distance from the CM.
Later, I'll calculate \bar R, the average <R> over time steps. Then I will get bending modulus via 
avg of <(R- \bar R)^2> = avg of <R^2> -2<R>\bar R + (\bar R)^2
*/
//for in-plane fluctuation calculation, need to subtract out global rotation and out-of-plane displacement

/*
double theta, phi;
double avg_theta = 0.;
double avg_phi = 0.;
*/
double distance_from_cm[DIMENSION];
double radial_projection, current_radius, current_radius2, in_plane_fluc_mag, inv_current_radius;
double total_current_radius = 0.;
double total_current_radius2 = 0.;
double inv_num_monos = 1./(double)(num_shell_monos-3);//subtracting out the 3 fixed monos
double in_plane_fluc[DIMENSION];
double total_in_plane_fluc_mag = 0.;
//calculate_cm();// getting ridiculous not to have a subroutine
//actually, not needed b/c I recalculate shell_cm in update_system

//first calculate total rotation that has occurred
//nope

//int wait;

//bending modulus calc
for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
if(mono_list[ii].get_update_bool())
{
 current_radius2 = 0.;
 radial_projection = 0.;
 in_plane_fluc_mag = 0.;
 for(kk = 0; kk < DIMENSION; kk++)
 {
	distance_from_cm[kk] = mono_list[ii].get_prev_pos(kk) - shell_cm[kk];
	current_radius2 += distance_from_cm[kk]*distance_from_cm[kk];
        in_plane_fluc[kk] = distance_from_cm[kk] - mono_list[ii].get_xbar(kk);
	radial_projection += in_plane_fluc[kk] * distance_from_cm[kk];
	in_plane_fluc_mag += in_plane_fluc[kk]*in_plane_fluc[kk];
 }

//fprintf(stderr, "dcm %g %g %g r2 %g inplane %g %g %g radial %g inmag %g invnum %g\n", distance_from_cm[0], distance_from_cm[1], distance_from_cm[2], current_radius2, in_plane_fluc[0], in_plane_fluc[1], in_plane_fluc[2], radial_projection, in_plane_fluc_mag, inv_num_monos);

 current_radius = sqrt(current_radius2);
 inv_current_radius = 1. / current_radius;
 in_plane_fluc_mag = sqrt(in_plane_fluc_mag);
 total_current_radius += current_radius;
 total_current_radius2 += current_radius2;
 radial_projection *= inv_current_radius / (current_radius * in_plane_fluc_mag);

 for(kk = 0; kk < DIMENSION; kk++)
 {
	in_plane_fluc[kk] -= radial_projection * distance_from_cm[kk];
	in_plane_fluc[kk] *= in_plane_fluc[kk];
        total_in_plane_fluc_mag += in_plane_fluc[kk]; 
 }//for kk


//fprintf(stderr, "2nd line dcm %g %g %g r2 %g inplane %g %g %g radial %g inmag %g invnum %g totr %g totr2 %g totalinmag %g\n", distance_from_cm[0], distance_from_cm[1], distance_from_cm[2], current_radius2, in_plane_fluc[0], in_plane_fluc[1], in_plane_fluc[2], radial_projection, in_plane_fluc_mag, inv_num_monos, total_current_radius, total_current_radius2, total_in_plane_fluc_mag);

//fprintf(stderr, "if update ii %i %g %g %g\n", ii, total_current_radius, total_current_radius2, total_in_plane_fluc_mag);

}//if update bool 

//fprintf(stderr, "end for loop ii %i %g %g %g\n", ii, total_current_radius, total_current_radius2, total_in_plane_fluc_mag);

}//for loop over shell monos

//fprintf(stderr, "before inv num %g %g %g\n", total_current_radius, total_current_radius2, total_in_plane_fluc_mag);
//scanf("%i", &wait);

//for bending mod calculation
total_current_radius *= inv_num_monos;
total_current_radius2 *= inv_num_monos;
total_in_plane_fluc_mag *= inv_num_monos;

//forget it, i'm going to fix a few monomers for this measurement to avoid rotations.
//perform rotation on x0

//fprintf(stderr, "%g %g %g\n", total_current_radius, total_current_radius2, total_in_plane_fluc_mag);
fprintf(fluctuation_file, "%llu %g %g %g\n", step, total_current_radius, total_current_radius2, total_in_plane_fluc_mag);
fflush(fluctuation_file);

}//end calculate_moduli_variables()


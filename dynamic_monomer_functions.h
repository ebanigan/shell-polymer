#include "dynamic_object_functions.h"


void dynamic_monomer::outward_pressure()
{
	double rad_vec[DIMENSION];
	int kk;
	double sum_sq = 0.;
	double norm;
	double motion[DIMENSION];
	for(kk = 0; kk < DIMENSION; kk++)
	{
		rad_vec[kk] =  get_prev_pos(kk) - shell_cm[kk];
		sum_sq = sum_sq + rad_vec[kk]*rad_vec[kk];
	}
	norm = 1./sqrt(sum_sq);
	for(kk = 0; kk < DIMENSION; kk++)
	{
		rad_vec[kk] = rad_vec[kk] * norm;
		motion[kk] = OS_PRESSURE*rad_vec[kk]*dt*get_tdiffusion_coeff()*invKT;
	}

	move(motion);
}




void dynamic_monomer::sys_force(dynamic_monomer *mono)
{
     double force;
     double separation2;
     double movement1[DIMENSION], movement2[DIMENSION];
     separation2 =  calculate_distance_sq(mono, PREVIOUS);
     int kk;
     for(kk = 0; kk < DIMENSION; kk++)
     {
/*if(distance condition) statements for now can remain inside the harmonic_repulsion function
if I add multiple interaction potentials, it may be good to calculate the distance in sys_force
and feed it into each interaction function*/
      force = harmonic_repulsion(mono, separation2, kk);

      movement1[kk] = (get_tdiffusion_coeff())*dt*force;//I used to multiply by invKT here, but that canceled out the temperature dependence-mainly because when I started playing with temperature, I took the temperature dependence out of the stiffness coefficients.160405
      movement2[kk] = -1.*dt*force*(*mono).get_tdiffusion_coeff();
     }

     move(movement1);
     (*mono).move(movement2);

     //fprintf(stderr, "mvmnt %g %g %g, SPRINGCONST %g\n", movement1[0], movement1[1], movement1[2], SPRINGCONST);

}//end of sys_force()



double dynamic_monomer::harmonic_repulsion(dynamic_monomer *mono, double separation2, int dd)
{
       double repulsion = 0.;
       double separation;
       int kk; 
       double tot_overlap_dist = (overlap_dist + (*mono).get_overlap_dist());
       double tot_overlap_dist2 = tot_overlap_dist*tot_overlap_dist;
       double effective_spring_const;
 
      if(separation2 < tot_overlap_dist2)
      {
        separation = sqrt(separation2);
        //overlap = sqrt(separation2) - diameter;
//really, we are compressing to springs (in series).  to get an effective spring constant to be SPRINGCONST (or whatever fed-in variable) as is conventional in my codes, need to multiply by factor of 2 in front 
	effective_spring_const = 2.*exc_vol_spring_const*(*mono).get_exc_vol_spring_const() / (exc_vol_spring_const + (*mono).get_exc_vol_spring_const());

        repulsion -= effective_spring_const*(calculate_1d_sep(mono, dd, PREVIOUS))*(1.-tot_overlap_dist/separation);

        //in general, springconst could be a property of each disk...and remember we're assuming all disks have same radius
        /*If calculate_1d_sep() > 0, RHS < 0 (since separation < diameter). So we want to subtract.
          In other words, (for ex.) x position of this object is greater than x position of mono,
          we want a positive repulsion, which requires a negative sign in front of the whole RHS (via subtraction).*/ 
        //note, other coefficients are already accounted for in sys_force() function to which this returns
      }

//fprintf(stderr, "repulsion = %g???\n", repulsion); 
  return repulsion; 
}//end harmonic repulsion


void dynamic_monomer::sys_torque(dynamic_monomer *mono)
{
     double torque= 0.;
     int kk;
     for(kk = 0; kk < DIMENSION; kk++)
     {
      if(calculate_distance_sq(mono, PREVIOUS) < diam2)
       torque = 0.;

      polarization[kk] += get_rdiffusion_coeff()*dt*invKT*torque;//?
     }     
}//end of sys_torque()


bool dynamic_monomer::proximity(dynamic_monomer *mono2)
{
    double test_dist2 = calculate_distance_sq(mono2, PREVIOUS);
    if(MONO_DIAM2 < test_dist2)
       return false;
    else if(MIN_DIST_FOR_POLYM_SQ > test_dist2)
       return false;
    else
       return true;
}//end of proximity()


#if (DIMENSION != 1)
/*NOTE: This condition will be the same for monomers added via barbed end polymerization
 as for monomers added to the side of a filament via branching.  This is because we ALWAYS
 want the incoming monomer to have a polarization vector (nearly) equal to the new bond 
 vector!  We only want a 70 degree angle between the "branch base" bond and the filament
 (noting that the polarizations of monomers in the main filament align with the filament
 bonds)*/

bool dynamic_monomer::orientation(dynamic_monomer *mono2)
{
     double cos_angle = 0.;//cosine of the angle between axes of orientation
     double bond_vector[DIMENSION];
     double bond_vector_mag = 0.;
     int kk;
     
     for(kk = 0; kk < DIMENSION; kk++)
     {
       bond_vector[kk] = -1.*calculate_1d_sep(mono2, kk, PREVIOUS);
       cos_angle += ((*mono2).get_prev_polarization(kk)) * bond_vector[kk];
       bond_vector_mag += (bond_vector[kk])*(bond_vector[kk]);
     }
     
     bond_vector_mag = sqrt(bond_vector_mag);
     cos_angle *= 1./bond_vector_mag;
    
     if(fabs(cos_angle - COS_F_ANGLE) > ORIENTATION_CUTOFF)
      return false;
     else
      return true;
}//end of orientation()


bool dynamic_monomer::incoming_angle(dynamic_monomer *incoming_mono, dynamic_monomer *minus_mono, bool for_a_branch)
{
     double test_angle = calculate_cosbond_angle(minus_mono, incoming_mono);
     //calculate_cosbond_angle calculates cosine of the bond angle!

//need to fix this so branch and filament nucleations have same angle range...
     if(!for_a_branch)
     {
       if(fabs(test_angle - COS_F_ANGLE) > INCOMING_ANGLE_CUTOFF)
         return false;
       else
         return true;
     }
     else
     {
       if(fabs(test_angle - COS_B_ANGLE) > INCOMING_ANGLE_CUTOFF)
         return false;
       else
         return true;
     }
}//end of incoming_angle()
#endif       


double dynamic_monomer::calculate_cosbond_angle(dynamic_monomer *minusptr, dynamic_monomer *second_mono)
{     
      int kk;
      double first_bond_vector[DIMENSION], second_bond_vector[DIMENSION];
      double first_vector_mag = 0.;
      double second_vector_mag = 0.;
      double dot_product = 0.;
      double costheta;//cosine of angle between first bond and second bond, to be calculated via dot product

      for(kk = 0; kk < DIMENSION; kk++)
      { 
        first_bond_vector[kk] = calculate_1d_sep(minusptr, kk, PREVIOUS);
        second_bond_vector[kk] = -1.*calculate_1d_sep(second_mono, kk, PREVIOUS);      
        
        first_vector_mag += (first_bond_vector[kk])*(first_bond_vector[kk]);
        second_vector_mag += (second_bond_vector[kk])*(second_bond_vector[kk]);
            
        dot_product += (first_bond_vector[kk])*(second_bond_vector[kk]);
      }
         
      first_vector_mag = sqrt(first_vector_mag);
      second_vector_mag = sqrt(second_vector_mag);
        
      costheta = dot_product / first_vector_mag / second_vector_mag;
      
      return costheta;
}//end of calculate_cosbond_angle()


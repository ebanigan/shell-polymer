#include "dynamic_monomer_functions.h"
//#include "neighbor_list.h"

//using namespace std;

//#include <vector>

/*
monomer_pair is an object composed of two pointers to monomer. These monomers
are connected to each other in a polymer/network.  In this way, the simulation
keeps track of every pairwise polymer/network bond interaction.
*/

class monomer_pair{
  protected:
//    long pair_id;
/*Relative position variables.  These can be calculated once per time step with
the update_relative_positions function and then used by each of the interaction 
functions.  Since these will be used by every pair on every time step, it is 
more efficient than calculating each as needed (which would be more than once).*/
      double separation, separation2, inv_separation;
      double relative_position[DIMENSION]; //position1 - position2
      double bond_strength, bend_strength;
      double orient_strength, nn_orient_strength;//orient_strength corresponds to "POLARIZABILITY," the strength of the interaction between polarization vectors and bond vectors.  nn_orient_strength corresponds to "POLARIZATION_INT," the strength of the interaction between nearby polarization vectors.
  
	double bond_length;  
  public:    
    //ids of monomers in pair
    dynamic_monomer *first;//first monomer's plus pointer points to second monomer
    dynamic_monomer *second;//second monomer's minus points to first monomer
    
    bool nearest_neighbors;//needed to determine preferred angle in polarization_interaction()
    bool polymer_pair;//indicates whether pair is for the interior polymer or not

    //monomer_pair(){}
    //constructor
    monomer_pair(dynamic_monomer *mono1, dynamic_monomer *mono2, bool polymer, double bond, double bend, double polarize, double nn_pol)
    {
       int kk;
       
       first = mono1;
       second = mono2;
       //set nearest_neighbors bool
       if((*second).minus_points_to() == (*first).get_id())
         nearest_neighbors = true;
       else
         nearest_neighbors = false;
       //interaction strength -- may be different, e.g., when comparing a crosslink pair to a ParA filament pair
       bond_strength = bond;
       bend_strength = bend;
       orient_strength = polarize;
       nn_orient_strength = nn_pol;
       polymer_pair = polymer;

	update_relative_positions();

	if(SPRINGS_ONLY || VARIABLE_BOND_LENGTH)
	  bond_length = sqrt((*second).calculate_distance_sq(first, PREVIOUS));
	else
	  bond_length = MONO_DIAM;

//fprintf(stderr, "bl %g ids %i %i\n", bond_length, (*first).get_id(), (*second).get_id());
    }
    
    //to be called at the end of update_system()
    void update_relative_positions()
    {
       int kk;
       separation2=0.;
       for(kk = 0; kk < DIMENSION; kk++)
       {
          relative_position[kk] = (*first).calculate_1d_sep(second, kk, PREVIOUS);
	  separation2 += relative_position[kk]*relative_position[kk];
       }
       separation = sqrt(separation2);
       inv_separation = 1./separation;
    }

    bool get_polymer_pair(){return polymer_pair;}
    double get_bond_length(){return bond_length;}
    void set_bond_length(double input){bond_length = input;}
    double get_bond_strength(){return bond_strength;}
    void set_bond_strength(double input){bond_strength = input;}
    double get_bend_strength(){return bend_strength;}
    double get_orient_strength(){return orient_strength;}
    double get_nn_orient_strength(){return nn_orient_strength;}

    double get_relative_position(int kk){return relative_position[kk];}
    double get_separation2(){return separation2;}
    double get_separation(){return separation;}
    double get_inv_separation(){return inv_separation;}
    
    
/*Function for the bond attraction between monomers.  Strength depends on distance between monomers.*/
    void bond_potential()
    {
     double attraction1[DIMENSION], attraction2[DIMENSION];//attractive bond force for "first" and "second", respectively
     double prefactor1, prefactor2;
     int kk;
     
     for(kk = 0; kk < DIMENSION; kk++)
     {
       attraction1[kk] = 0.;
       attraction2[kk] = 0.;
     }

// the idea is that I want shell springs to be repulsive (in addition to the exc. vol of nodes)
     if((separation2 > bond_length*bond_length) || ((SPRINGS_ONLY && !EXTENSIONAL_SPRINGS_ONLY) || (VARIABLE_BOND_LENGTH && ((*first).get_id() >= NUMBER_IN_POLYMER))) )//if variable shell bond length, we need to make sure polymer mono pairs are treated properly (i.e., no repulsion here)
//since bond_length will generally be less than MONO_DIAM, this potential isn't the main source of repulsive interactions
     {
      //this expression may have to be changed if we add hydrodynamics, since tdiffusion coefficient will change
      prefactor1 = ((*first).get_tdiffusion_coeff())*dt*bond_strength;//no invKT here! -- unless stiffness coeffs get their KT back.. 160405
      prefactor2 = ((*second).get_tdiffusion_coeff())*dt*bond_strength;
      /*if (bond_strength < 99.) 
         fprintf(stderr, "bond %g, pref %g\n", bond_strength,prefactor1); 
      else if (bond_strength > 101.)
         fprintf(stderr, "bond %g, pref %g\n", bond_strength, prefactor1);*/
      for(kk = 0; kk < DIMENSION; kk++)
      {     
       attraction1[kk] -= prefactor1*(relative_position[kk])*(1.-bond_length*inv_separation);
       //since separation > diameter, when calc_1d_sep > 0, RHS > 0... so subtract
       //e.g., if first's x is greater than second's x, we want first's x to decrease.
       attraction2[kk] += prefactor2*(relative_position[kk])*(1-bond_length*inv_separation);      
      }

//fprintf(stderr, "bond strength = %g length = %g\n", bond_strength, bond_length);
 
     (*first).move(attraction1);
     (*second).move(attraction2);
     }
    }//end of bond_potential(), which is generic for all pairs


#if (DIMENSION != 1)
/*Function for moving monomers based on angle between bond vectors (bending).  
This bend potential will work in 2D and 3D.  There could be different theta0's 
for polymer bond interactions if we want.
interaction. For more information, see D.C. Rapaport (2004) p. 284 -- "bond angle" 
potential. */
    void bending_potential(monomer_pair *pair_ij)
    {
         int ll;       
         double cos_angle;
         double rkj[DIMENSION], rji[DIMENSION];
         double rkjmag2= separation2;
         double rjimag2= (*pair_ij).get_separation2();
         double rji_times_rkj;
         double rji_dot_rkj = 0.;
         double inv_rji2rkj2;
         double prefactor1, prefactor2, prefactor3;
         double bend_move_i[DIMENSION], bend_move_j[DIMENSION], bend_move_k[DIMENSION];

         for(ll = 0; ll < DIMENSION; ll++)
         {
            rkj[ll] = -1.*relative_position[ll];
            rji[ll] = -1.*(*pair_ij).get_relative_position(ll);
            rji_dot_rkj += (rkj[ll])*(rji[ll]);
         }
         inv_rji2rkj2 = 1./(rkjmag2*rjimag2);
         rji_times_rkj = (*pair_ij).get_separation()*separation;
         
         cos_angle = rji_dot_rkj/rji_times_rkj;

           prefactor1 = (*first).get_tdiffusion_coeff()*dt*bend_strength*(COS_F_ANGLE - cos_angle);
           prefactor2 = (*second).get_tdiffusion_coeff()*dt*bend_strength*(COS_F_ANGLE - cos_angle);
           prefactor3 = (*((*pair_ij).first)).get_tdiffusion_coeff()*dt*bend_strength*(COS_F_ANGLE - cos_angle);

         //for details about the bending force, see analytics notes from 7/4/08
         for(ll = 0; ll < DIMENSION; ll++)
         {
           bend_move_i[ll] = -1.*(rkj[ll])*rji_times_rkj + cos_angle*(rji[ll])*rkjmag2;
           bend_move_i[ll] *= inv_rji2rkj2;
           bend_move_i[ll] *= prefactor3;
  
           bend_move_k[ll] = (rji[ll])*rji_times_rkj - cos_angle*(rkj[ll])*rjimag2;
           bend_move_k[ll] *= inv_rji2rkj2;
           bend_move_k[ll] *= prefactor2;

           bend_move_j[ll] = (rkj[ll] - rji[ll])*rji_times_rkj - cos_angle*(rkjmag2*rji[ll]-rjimag2*rkj[ll]);
           bend_move_j[ll] *= inv_rji2rkj2;
           bend_move_j[ll] *= prefactor1;
         } 

         (*((*pair_ij).first)).move(bend_move_i);
         (*first).move(bend_move_j);
         (*second).move(bend_move_k);
    }//end of bending_potential()


/*Rotates monomers according to angle between polarization vector and (1 or 2) nearest bond vectors.
Uses potential = .5*K*(costheta - costheta0)^2.
*/
    void orientation_potential()
    {
         double bond_vector[DIMENSION], bhat[DIMENSION];
//         double bond_vector_mag = separation;
         double inv_bond_vector_mag = inv_separation;
//         double inv_bond_vector_mag2;
         double pol1[DIMENSION], pol2[DIMENSION];
         double rot1[DIMENSION], rot2[DIMENSION];
         double move1[DIMENSION], move2[DIMENSION];
         double costheta0_1, costheta0_2;
         double p1_dot_bhat = 0.;
         double p2_dot_bhat = 0.;
         double r_prefactor1 = orient_strength*((*first).get_rdiffusion_coeff())*dt;
         double r_prefactor2 = orient_strength*((*second).get_rdiffusion_coeff())*dt;
         double t_prefactor1 = orient_strength*((*first).get_tdiffusion_coeff())*dt;
         double t_prefactor2 = orient_strength*((*second).get_tdiffusion_coeff())*dt;//not used below
         double temp1, temp2;
         
         //polarizations line up with bond vectors. 
          costheta0_1 = COS_F_ANGLE;
          costheta0_2 = costheta0_1;

         int kk;
         
         for(kk = 0; kk < DIMENSION; kk++)
         {
           bond_vector[kk] = -1.*relative_position[kk];
           bhat[kk] = bond_vector[kk]*inv_bond_vector_mag;
           pol1[kk] = (*first).get_prev_polarization(kk);
           pol2[kk] = (*second).get_prev_polarization(kk);
           p1_dot_bhat += (pol1[kk])*bhat[kk];
           p2_dot_bhat += (pol2[kk])*bhat[kk];
           rot1[kk] = 0.;
           rot2[kk] = 0.;
           move1[kk] = 0.;
           move2[kk] = 0.;
         }
         
         //Calculation of rotation of first and second monmer polarizations
         for(kk = 0; kk < DIMENSION; kk++)
         {
//minus sign from -del U = F incorporated by reversing order of costheta0 and p_dot_bhat
           rot1[kk] = r_prefactor1*(costheta0_1-p1_dot_bhat)*(bhat[kk]-p1_dot_bhat*pol1[kk]);
           rot2[kk] = r_prefactor2*(costheta0_2-p2_dot_bhat)*(bhat[kk]-p2_dot_bhat*pol2[kk]);
         }

         //feedback into bond
         for(kk = 0; kk < DIMENSION; kk++)
         {
           //effect on bond due to polarization 1
           temp1 = inv_bond_vector_mag*(pol1[kk]-p1_dot_bhat*bhat[kk])*(costheta0_1 - p1_dot_bhat);
           move2[kk] += t_prefactor2*temp1;
           move1[kk] -= t_prefactor1*temp1;

           //effect on bond due to polarization 2
           temp2 = inv_bond_vector_mag*(pol2[kk]-p2_dot_bhat*bhat[kk])*(costheta0_2 - p2_dot_bhat);
           move2[kk] += t_prefactor2*temp2;
           move1[kk] -= t_prefactor1*temp2;
         }


         (*first).move(move1);
         (*second).move(move2);
         (*first).rotate(rot1);
         (*second).rotate(rot2);
    }
#endif
#if CROSSLINK
void crosslink_potential()
{
     //attractive (or repulsive) force due to crosslinker attraction for "first" and "second", respectively
	 double motion1[DIMENSION], motion2[DIMENSION];
	 double prefactor1, prefactor2;
         int kk;
     
     for(kk = 0; kk < DIMENSION; kk++)
     {
       motion1[kk] = 0.;
       motion2[kk] = 0.;
     }
     
	 #if CROSSLINK_ATTRACTIVE_ONLY
     if(separation2 > CROSSLINK_LENGTH2)
     {
	 #endif
      //this expression may have to be changed if we add hydrodynamics, since tdiffusion coefficient will change
      prefactor1 = ((*first).get_tdiffusion_coeff())*dt*CROSSLINK_SPRING;
      prefactor2 = ((*second).get_tdiffusion_coeff())*dt*CROSSLINK_SPRING;
      for(kk = 0; kk < DIMENSION; kk++)
      {     
       motion1[kk] -= prefactor1*(relative_position[kk])*(1.-CROSSLINK_LENGTH*inv_separation);
       //since separation > diameter, when calc_1d_sep > 0, RHS > 0... so subtract
       //e.g., if first's x is greater than second's x, we want first's x to decrease.
       motion2[kk] += prefactor2*(relative_position[kk])*(1.-CROSSLINK_LENGTH*inv_separation);
      }
	 #if CROSSLINK_ATTRACTIVE_ONLY
     }
	 #endif

     (*first).move(motion1);
     (*second).move(motion2);
}
#endif



/*
 * Define RANDOM_SIDE_BIND_SITES to determine whether or not there is any ordering to the locations of side binding sites.
 * If SIDE_BIND_RING is false, then these will be called in filament_interactions() if RANDOM_SIDE_BIND_SITES is false.
 *
 * This potential, coupled with the orientation_potential() called when SIDE_BINDING_RING is false, are sufficient to identify a 
 * distinct, preferred location on the monomer for a binding site.  When this potential is absent, while there will still be a 
 * single binding site, given by the polarization vector, the site will diffuse freely in the ring perpendicular to the bond vector.
 * */

/*
 * Added the bool nearest_neighbors to monomer_pair.  Also added vector<monomer_pair> nextnnpairs to code.
 * Then, in filament_interactions(), call next_nearest_neighbor_list.polarization_interaction(), as well as pairs.polarization_interaction().
 * */

//potential that tends to align polarization vectors at POL_ANGLE
void polarization_interaction()
{
  double costheta0;
  double interaction_strength = nn_orient_strength;
  if(nearest_neighbors)
  {
    costheta0 = COS_POL_ANGLE;//nearest neighbors want to have POL_ANGLE between polarizations
  }//end of if(nearest neighbors)
  else
  {
    costheta0 = COS_TWO_POL_ANGLE;//next nearest neighbors prefer 2*POL_ANGLE between polarizations
  }

         double pol1[DIMENSION], pol2[DIMENSION];
         double rot1[DIMENSION], rot2[DIMENSION];
         double p1_dot_p2 = 0.;
         double r_prefactor1 = interaction_strength*((*first).get_rdiffusion_coeff())*dt;
         double r_prefactor2 = interaction_strength*((*second).get_rdiffusion_coeff())*dt;

         int kk;
         
         for(kk = 0; kk < DIMENSION; kk++)
         {
           pol1[kk] = (*first).get_prev_polarization(kk);
           pol2[kk] = (*second).get_prev_polarization(kk);
           p1_dot_p2 += (pol1[kk])*(pol2[kk]);
           rot1[kk] = 0.;
           rot2[kk] = 0.;
         }
         
         //Calculation of rotation of first and second monmer polarizations
         for(kk = 0; kk < DIMENSION; kk++)
         {
//minus sign from -del U = F incorporated by reversing order of costheta0 and p1_dot_p2
           rot1[kk] = r_prefactor1*(costheta0-p1_dot_p2)*(pol2[kk]-p1_dot_p2*pol1[kk]);
           rot2[kk] = r_prefactor2*(costheta0-p1_dot_p2)*(pol1[kk]-p1_dot_p2*pol2[kk]);
         }

         (*first).rotate(rot1);
         (*second).rotate(rot2);
}//end of polarization_interaction()

};


extern vector<monomer_pair> pairs;//polymer bonds in backbone of polymer
extern vector<monomer_pair> crosslinkpairs;




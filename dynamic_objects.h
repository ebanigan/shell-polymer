#include "monomer_events.h"
//this file has the classes dynamic_object, dynamic_monomer


class dynamic_object{
      protected:
        double pos[DIMENSION];//array of object's position (1 coordinate per dimension)
        double prev_pos[DIMENSION];//position at end of previous timestep
        //double vel[DIMENSION];//velocity. not used (6/10/08)
        
        double polarization[DIMENSION];
        double prev_polarization[DIMENSION];
    
        //object properties
        //double friction;
        double radius, radius2, diameter, diam2, invdiam2;//monomer radius, diameter squared, 1/diam sq
        double mass;
        double tdiffusion_coeff, rdiffusion_coeff;
        double overlap_dist;
	double tdiff_stdev;

        double z_brownian_move;
//holds the value of the Brownian move in the z-direction so that we know how much of movement is due to actual forces at each step

      public:
      //sample constructor
      /*dynamic_object()
      {
        int kk;
        double sum_sq = 0.;
        
        for(kk = 0; kk < DIMENSION; kk++)
        {
          pos[kk] = unif_rand()*L[kk];
          prev_pos[kk] = pos[kk]; //just to give it an initial value
        }
                        
       for(kk = 0; kk < DIMENSION; kk++)
       {
         //want to uniformly select polarization vector, by picking point in sphere and choosing
         //orientation to be in direction of vector to that pt.
         polarization[kk] = 2*(unif_rand()-0.5);
         prev_polarization[kk] = polarization[kk];
         sum_sq += (polarization[kk])*(polarization[kk]);
       
         if(sum_sq > 1.)
         {
           sum_sq = 0.;
           kk = -1;
         }//if sum_sq > 1, outside sphere, try again
       }
       
       //need to normalize polarization[]
       for(kk = 0; kk < DIMENSION; kk++)
       {
          polarization[kk] *= 1./sqrt(sum_sq);
       }       

        radius = MONO_RAD;
        radius2 = radius*radius;
        overlap_dist = radius;
        tdiffusion_coeff = TDIFF_COEFF;
        rdiffusion_coeff = RDIFF_COEFF;
        diameter = 2.*radius;
        friction = FRICTION_COEFF;
        diam2 = diameter*diameter;// (2r)^2
        invdiam2 = 1./diam2;
        mass = MONO_MASS; 
        }*/
        
        virtual void rotational_brownian_move();
        virtual void translational_brownian_move();
        //virtual void sys_force(void) =0;
             
//functions return positions and orientations
        double get_pos(int dd){return pos[dd];}
        double get_prev_pos(int dd){return prev_pos[dd];}
        double get_polarization(int dd){return polarization[dd];}
        double get_prev_polarization(int dd){return prev_polarization[dd];}

//functions return properties of the monomer
        double get_diameter(){return diameter;}
        double get_radius(){return radius;}
        double get_radius2(){return radius2;}
        double get_tdiffusion_coeff(){return tdiffusion_coeff;}//returns the translational diffusion coefficient
                                  //this is a function because it needs to be calculated for hydrodynamics sims
                                  //end b/c we have different dynamic objs
        double get_rdiffusion_coeff(){return rdiffusion_coeff;}
	void set_tdiffusion_coeff(double input){tdiffusion_coeff = input; tdiff_stdev = sqrt(2.*tdiffusion_coeff*dt);}
	void set_tdiff_stdev(double input){tdiff_stdev = input;}//to be used only with low temperature initialization 
        void set_rdiffusion_coeff(double input){rdiffusion_coeff = input;}
        double get_overlap_dist(){return overlap_dist;}
	void set_radius(double input){radius = input; radius2 = radius*radius; diameter = 2.*radius; diam2 = diameter*diameter; invdiam2 = 1./diam2; overlap_dist = radius;}

//set prev_pos and _polarization set prev = current
        void set_prev_pos();//sets prev_pos = pos, to be used at the beginning of each time step
        void set_prev_polarization();//sets prev_polarization = polarization, to be used at beginning of time step

	double get_z_brownian_move(){return z_brownian_move;}

//set (current) pos and polarization set pos and polarization to an input
        void set_pos(int dd, double position_dd){pos[dd] = position_dd;}//sets position equal to some input
        void set_polarization(int dd, double directionmag){polarization[dd] = directionmag;}
        
        //returns angular orientation in polar coordinates
        double get_theta();//3d and 2d -- note that theta follows convention for 3d (azimuthal) and 2d meanings (polar)
        double get_phi();//3d only polar angle
        

//functions (nontrivial)
        //calculate velocity via 2pt. differentiation
      //double calculate_speed();


//should feed in two generic dynamical objects for these two:

//calculate distance squared between object and obj2 (fed into function), return distance2 (calculated)
//uses positions from beginning of time step (prev. pos)
      double calculate_distance_sq(dynamic_object *obj2, bool previous_positions);
//calculates difference in x, y, or z positons of 2 objs
      double calculate_1d_sep(dynamic_object *obj2, int dd, bool previous_positions);
      
      
      /*Keeps objs within boundaries of periodically bounded box.*/        
        void check_translational_pbc();
/*resets pol vect mag to 1 (maintains vector direction)*/
        void renorm_polarization();
        
        //moves obj by a specified amount
//could generalize translational brownian_move to just be move with an input of DIMENSION#-element array of gaussian1()
        void move(double movement[DIMENSION]);
	void move_1d(int dd, double mvmnt){pos[dd] += mvmnt;}

        void rotate(double rotation[DIMENSION]);



};//end of dynamic_object




/*dynamic_monomer is a monomer that moves and rotates (has dynamics).*/
class dynamic_monomer:public monomer, public dynamic_object{
      protected:
    
	bool load_mono;
	double sign_load;
	bool binding_mono;

	double load_rest_length;

	double exc_vol_spring_const;

        double x0[DIMENSION]; //intial position of monomer
	double xbar[DIMENSION]; // time avged position        

        int pair_num2, crosslink_num, nextnn_pair_num2;
/*
 * pair_num2 is id of "monomer_pair" in which monomer is  "second"
 * crosslink_num is the crosslink pair to which the monomer belongs
 * nextnn_pair_num2 is the id of the next nearest neighbor monomer_pair in which monomer is second.  applicable when there is an interaction between polarization vectors
 */
		

        bool update_bool;//tells sim whether or not to update this mono

 
      public:
        void outward_pressure();

        void set_update_bool(bool update_status){update_bool = update_status;}
        bool get_update_bool(){return update_bool;}


        //functions for getting and changing "pair_num' variables
        long get_pair_num2(){return pair_num2;}
        long get_crosslink_num(){return crosslink_num;}
        long get_nextnn_pair_num2(){return nextnn_pair_num2;}
        void set_pair_num2(long number){pair_num2 = number;}
	void set_crosslink_num(long number){crosslink_num = number;}
        void set_nextnn_pair_num2(long number){nextnn_pair_num2 = number;}
        void decrement_pair_num2(){pair_num2--;}
	void decrement_crosslink_num(){crosslink_num--;}
        void decrement_nextnn_pair_num2(){nextnn_pair_num2--;}

		

//default constructor
dynamic_monomer(){}

//useful constructor      
dynamic_monomer(long num, double invdrag)
{
        int kk;
        long jj;
        double sum_sq = 0.;


        exc_vol_spring_const = SPRINGCONST;       
        crosslinked = false; 
        shell_poly_link = false;
        id = num;
        update_bool = true;
        pair_num2 = -1;
	crosslink_num = -1;
        nextnn_pair_num2 = -1;
        depolymerize_prob = 0.;
        
        plus = NULL;
        minus = NULL;
        plus_is_polymerized = false;
        minus_is_polymerized = false;
		
	load_mono = false;
	sign_load = 0.; 
 	binding_mono = false;

        for(kk = 0; kk < DIMENSION; kk++)
        {
          set_pos(kk, unif_rand()*L[kk]);
          x0[kk] = pos[kk];
	  xbar[kk] = x0[kk];
        }
        set_prev_pos();//just to give it an initial value
                        
       for(kk = 0; kk < DIMENSION; kk++)
       {
         //want to uniformly select polarization vector, by picking point in sphere and choosing
         //orientation to be in direction of vector to that pt.
         polarization[kk] = 2*(unif_rand()-0.5);
         sum_sq += (polarization[kk])*(polarization[kk]);
       
         if(sum_sq > 1.)
         {
           sum_sq = 0.;
           kk = -1;
         }//if sum_sq > 1, outside sphere, try again
       }
       
       //need to normalize polarization[]
       for(kk = 0; kk < DIMENSION; kk++)
       {
          polarization[kk] *= 1./sqrt(sum_sq);
       }
       set_prev_polarization();

        radius = MONO_RAD;
        radius2 = radius*radius;
        overlap_dist = radius;
        tdiffusion_coeff = invdrag;
	tdiff_stdev = sqrt(2.*tdiffusion_coeff*dt);
        rdiffusion_coeff = RDIFF_COEFF/TDIFF_COEFF*invdrag;//drag = TDIFF/drag factor
        diameter = 2.*radius;
        //friction = KT/MONO_MASS/tdiffusion_coeff;
        diam2 = diameter*diameter;// (2r)^2
        invdiam2 = 1./diam2;
        mass = MONO_MASS;

}//end of constructor

double get_x0(int dd){return x0[dd];}    
void set_x0(){for(int kk = 0; kk < DIMENSION; kk++) x0[kk] = pos[kk];}

double get_xbar(int dd){return xbar[dd];}
void set_xbar(int dd, double xvec){xbar[dd] = xvec;}

double get_exc_vol_spring_const(){return exc_vol_spring_const;}
void set_exc_vol_spring_const(double input){exc_vol_spring_const = input;}
void choose_load_mono(double input, double dist_from_center)
{
  load_mono = true; 
  sign_load = input;
  if(LENGTH_CONTROLLED_LOAD)
    load_rest_length = dist_from_center - fabs(prev_pos[0] - 0.5*L[0]);
}
bool get_load_mono(){return load_mono;}
double get_sign_load(){return sign_load;}
bool get_binding_mono(){return binding_mono;}
void choose_binding_mono(){binding_mono = true;}
void turn_off_binding_mono(){binding_mono = false;}

//subroutine to pull apart shell
       void load()
       {
	 if(!LENGTH_CONTROLLED_LOAD)	
            pos[0] += get_sign_load()*get_tdiffusion_coeff()*dt*F_LOAD;
	 else
	 {
	     pos[0] += get_tdiffusion_coeff()*dt* PIPETTE_STIFFNESS*(0.5*L[0] + get_sign_load()*(pipette_dist_from_center - load_rest_length) - prev_pos[0]);
	 }
       }

double get_load_rest_length(){return load_rest_length;}


//adds stochastic motion to particle -- adds gaussian1() to each dimension of position
        void translational_brownian_move(unsigned long long step)
        {
          int kk;

          for(kk = 0; kk < DIMENSION; kk++)
          {
              pos[kk] += tdiff_stdev * gaussian_inverf(step);//gaussian_rand();//gaussian_std();
          }
        }//may want some monos to have higher drag. use "<" instead of "==" because these are doubles not longs... don't want weird numerical errors.


/*Rotates the monomer randomly with angle from gaussian distribution. In 3d, rotation is about axis chosen
from a uniform distribution. In 2d, rotation is about center of monomer. Does nothing for 1d monomers (6/20/08).
     method:
     1.choose rotation
     2.translate to alpha, phi
     3.apply rotation
     4.renormalize polarization vector (eliminate numerical error)*/
     //void rotational_brownian_move();


//adds systematic forces, for ex. a soft repulsion, adjusts position accordingly
        void sys_force(dynamic_monomer *mono);

        double harmonic_repulsion(dynamic_monomer *mono, double separation2, int dd);


//calculates sq displacement from initial position
      double calculate_sq_disp();
//calculates displacement from initial position
      double calculate_disp(int dd);
      

        //calculates angle between two bonds
        //no func to calc angle between orientation and bond
        double calculate_cosbond_angle(dynamic_monomer *minusptr, dynamic_monomer *second_mono);
        
};


extern vector<dynamic_monomer> mono_list;
extern dynamic_monomer central_mono;



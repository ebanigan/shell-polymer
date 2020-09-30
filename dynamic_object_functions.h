#include "dynamic_objects.h"


void dynamic_object::translational_brownian_move()
{    
    int kk;

    for(kk = 0; kk < DIMENSION; kk++)
    {
	pos[kk] += gaussian_std();
    }



//    check_translational_pbc();   
}//end of translational_brownian_move()


/*
void dynamic_object::rotational_brownian_move()
{

for(int kk = 0; kk < DIMENSION; kk++)
 polarization[kk] += gaussian2();

renorm_polarization();
}
*/


void dynamic_object::rotational_brownian_move()
{
/*method:
     1.choose rotation
     2.translate to alpha, phi
     3.apply rotation
     4.renormalize polarization vector (eliminate numerical error)*/
     
    int kk;   
    double nhat[DIMENSION];
    double nhatdotpolarization = 0.;
    double sum_sq = 0.;
    double alpha;
    double c_alpha, s_alpha;
     
    //choose number uniformly inside sphere in 3D, in lower D skip this step.
    #if (DIMENSION == 3)
     for(kk = 0; kk < DIMENSION; kk++)
     {

       nhat[kk] = 2.*(unif_rand()-0.5);

       sum_sq += (nhat[kk])*(nhat[kk]);
       
       if(sum_sq > 1.)
       {
         sum_sq = 0.;
         kk = -1;
       }
     }//choose polarization of rotation uniformly
 
     for(kk = 0; kk < DIMENSION; kk++)
       nhat[kk] *= 1./sqrt(sum_sq);
    
     for(kk = 0; kk < DIMENSION; kk++)
       nhatdotpolarization += (nhat[kk])*prev_polarization[kk];  
    #endif

if(rdiffusion_coeff < .99*RDIFF_COEFF)
     alpha = gaussian6();
else     
     alpha = gaussian2();//gaussian choice of alpha// PI*sqrt(sum_sq);//uniform choice of alpha
     c_alpha = cos(alpha);
     s_alpha = sin(alpha);
     
     #if (DIMENSION == 3)
      for(kk = 0; kk < DIMENSION; kk++)
      {
       polarization[kk] += (prev_polarization[kk])*(c_alpha-1)+(nhat[kk])*nhatdotpolarization*(1.-c_alpha)
                 + s_alpha * ((prev_polarization[(kk+1)%DIMENSION])*nhat[(kk+2)%DIMENSION]
                 - (prev_polarization[(kk+2)%DIMENSION])*nhat[(kk+1)%DIMENSION]);
      }/*From Wolfram Mathworld "Rotation Formula."  Also see Goldstein 1980, p164; Varshalovich et al., p.24.
       This step used to have somewhat large numerical error. Not sure if it still does.
       So, renormalize "polarization" every time (at the end of a time step).
       9/23/08 note that the scheme has been slightly modified (+= instead of =) so that we can use potentials
       to modify polarization[kk] in other functions. In particular, we add p'-p_prev to p_current (where p'
       is p_prev under the rotation.*/

     #elif (DIMENSION == 2)
       polarization[0] = c_alpha*prev_polarization[0] + s_alpha*prev_polarization[1];
       polarization[1] = c_alpha*prev_polarization[1] - s_alpha*prev_polarization[0];
     #endif//For serious 2D system, it would be easier to just keep track of angle 
      //from x polarization (theta) and shift that by Gaussian perturbations

//Taken out 9/23/08
//     renorm_polarization();
}//end of rotational_brownian_move()



void dynamic_object::set_prev_pos()
{
     for(int kk = 0; kk < DIMENSION; kk++)
       prev_pos[kk] = pos[kk];
}

void dynamic_object::set_prev_polarization()
{
     renorm_polarization();
     for(int kk = 0; kk < DIMENSION; kk++)
       prev_polarization[kk] = polarization[kk];
}


double dynamic_object::get_theta()
{
     double theta;
     
#if (DIMENSION == 3)
       theta = acos(polarization[2]);//no adjustments needed since 0<theta<pi in 3d
#elif (DIMENSION == 2)
       theta = acos(polarization[0]);
       if(polarization[1] < 0.)
         theta = TWOPI - theta;
#else
       theta = -999.;
#endif
     return theta;
}//end of get_theta()


double dynamic_object::get_phi()
{       
       double phi;
       phi = acos((polarization[0]) / sin(get_theta()));
       if(polarization[1] < 0.)
        phi = TWOPI - phi;
       
       return phi;
       /*     
     polarization0= stheta cphi
     polarization1= stheta sphi
     polarization2= ctheta*/
}//end of get_phi()

/*
double dynamic_object::calculate_speed()
{
   int kk;
   double speed = 0.;
   double temp;
   
   for(kk = 0; kk < DIMENSION; kk++)
   {
     temp = pos[kk] - prev_pos[kk];
     if(temp > (L[kk])/2.)
       temp -= L[kk];
     if(temp < -(L[kk])/2.)
       temp += L[kk];
       
     speed += temp*temp;
   }
     
   speed = invdt*sqrt(speed);
   
   return speed;
}//end calc speed
*/

double dynamic_object::calculate_distance_sq(dynamic_object *obj2, bool previous_positions)
{
      double temp;
      int kk;
      
      double distancesq = 0.;
      
      for(kk = 0; kk < DIMENSION; kk++)
      {

        temp = calculate_1d_sep(obj2, kk, previous_positions);//uses positions from beginning of time step

        distancesq += temp*temp;
      }
      return distancesq;
}//end of caluclate_distance_sq


double dynamic_object::calculate_1d_sep(dynamic_object *obj2, int dd, bool previous_positions)
{
     double separation;

     if(previous_positions)
      separation = prev_pos[dd] - (*obj2).get_prev_pos(dd);
     else
      separation = pos[dd] - (*obj2).get_pos(dd);

      if(separation < -0.5*L[dd])
          separation += L[dd]; 
      else if(separation > 0.5*L[dd])
          separation -= L[dd];

     return separation;
}//end of calculate_1d_sep


void dynamic_object::check_translational_pbc()
{
   int kk;        
   for(kk = 0; kk < DIMENSION; kk++)
   {
if(!FINITE_BOX)
{
     if(pos[kk] < 0.)
       pos[kk] += L[kk];
     else if(pos[kk] >= L[kk])     
       pos[kk] -= L[kk];
}
else
{
	if(pos[kk] < 0.)
	  pos[kk] = 0.00001;
	else if(pos[kk] > L[kk])
	  pos[kk] = 0.99999*L[kk];//numerical error issue?
}
/*To be perfectly accurate, this could be a while statement, but none of the particles 
should move that much in one time step!*/
   }
}//end of check_translational_pbc



void dynamic_object::renorm_polarization()
{
     int kk;
     double polsum_sq = 0.;
     for(kk = 0; kk < DIMENSION; kk++)
       polsum_sq += (polarization[kk])*polarization[kk];
     for(kk = 0; kk < DIMENSION; kk++)
       polarization[kk] *= 1./sqrt(polsum_sq);
}


void dynamic_object::move(double movement[DIMENSION])
{
    int kk;
    for(kk = 0; kk < DIMENSION; kk++)
      pos[kk] += movement[kk]; 

//    check_translational_pbc();
}//end move()



void dynamic_object::rotate(double rotation[DIMENSION])
{
    int kk;
    for(kk = 0; kk < DIMENSION; kk++)
      polarization[kk] += rotation[kk];
    
//    renorm_polarization();
}//end of rotate()




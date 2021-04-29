#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include "RNG_taus.h"

#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/detail/erf_inv.hpp>
//#include "../libraries/boost/libs/math/include/boost/math/special_functions/erf.hpp"
//#include "../libraries/boost/libs/math/include/boost/math/special_functions/erf_inv.hpp"
using namespace boost::math;


#define COMPRESS cl_compress
#define WALLSPRING 100.
#define WALLTHICKNESS 1.
#define WALLDIFF 0.001

//options relevant for making movies
#define LOWT_PROTOCOL cl_lowt_protocol
#define BOX_RESIZE cl_box_resize
#define EXTENDED_RELAXATION cl_extended_relaxation

#define CALC_MODULI cl_calc_moduli
#define HYSTERESIS cl_hysteresis

#define PRINT_YLM false
#define PRINT_SHELL_SPRING false//movie
#define PRINT_FOURIER true

#define LOAD_TEST false
#define LOAD_TEST_SPRING_CONST 15.

#define FINITE_BOX true

#define SPRINGS_ONLY cl_springs_only
#define EXTENSIONAL_SPRINGS_ONLY cl_extensional_springs_only
#define THERMAL cl_thermal
#define LENGTH_DEPENDENT_SPRINGS cl_length_dependent_springs
#define VARIABLE_BOND_LENGTH cl_variable_bond_length
#define LDEP_FACTOR cl_ldep_factor

#define ELLIPSOID cl_ellipsoid
#define MAJOR_AXIS (1.4*SHELL_RADIUS)
#define MINOR_AXIS_1 SHELL_RADIUS
#define MINOR_AXIS_2 MINOR_AXIS_1

#define CYLINDER cl_cylinder
#define CYLINDER_LENGTH cl_cylinder_length

#define MODIFIED_SHELL_EXC_VOL cl_modified_shell_exc_vol

#define TRIALNUMBER cl_trialnumber
#define rseed cl_rseed//419//random # seed
#define WRITE true
#define RESTART cl_restart

#define LOAD false


#define DNA_DRAG_FACTOR (1./TDIFF_COEFF)
#define DRAG_ALL_SITES false


#define SHELL_RADIUS cl_shell_radius//5.0
#define NUM_RELAXATION_STEPS 200000//((long)(50e3))
#define num_shell_monos cl_num_shell_monos//((int)(16.*SHELL_RADIUS*SHELL_RADIUS))//1.4
#define LOWER_BOUND_CONNECTIVITY cl_lower_bound_connectivity//3//4
#define UPPER_BOUND_CONNECTIVITY cl_upper_bound_connectivity//7

#define NUM_BINDING_MONOS ((int)(0.0*num_shell_monos))
#define NUM_LOAD_MONOS cl_num_load_monos
#define OS_PRESSURE cl_os_pressure
#define SHELL_BOND_SPRING cl_shell_bond_spring//(1.0*BOND_SPRING)
#define SOLID_INTERIOR cl_solid_interior

#define NUMBER_OF_MONOMERS (num_shell_monos+NUMBER_IN_POLYMER)//total number of "monomer" objects in the system
#define INITIALIZE_DIMERS true 
/*Set to true to initialize the system with dimers. 
If set to false, spontaneous nucleation is allowed*/

#define BOOK_KEEPING_TEST true//set to true to run book_keeping_tests() -- need this true to run number_free()
#define PRINT_NUM_IN_FIL true//set to true if want to print #monos in filaments (instead of free + or - ends) to file

#define HARMONIC_REPEL true/*True: Disks will repel each other with a soft harmonic repulsion*/
#define ROTATE false

#define ORPOT true//bool for orientatin potential on/off
#define SINGLE_POLYMER false //True: Initialize one long, straight chain of monomers
#define PDB cl_pdb//false//true

#define dt cl_dt//5e-4//2.5e-4 for high temp or l-dependent springs
#define NUMSTEPS cl_numsteps//((long)30000000+1)//((long)20e6+1)
#define NUMSKIP 1000//2000

#define NUM_INIT_STEPS 10//25000//((long)2e6)//10e6
#define NUM_INIT_POLY_STEPS (10000*NUM_INIT_STEPS)

#define DIMENSION 3
//length of sim box in each dimension... in general these should be greater than 3.3*MONO_DIAM so that NX,NY,NZ >= 3
#define LX cl_lx//50.//(200*MONO_DIAM)
#define LY cl_ly//(4.05*SHELL_RADIUS)//(200*MONO_DIAM)
#define LZ LY//(200*MONO_DIAM)

#define KT cl_kt//1.e-3//low T = 1.e-3 //high T = 10.//1.
//#define VISCOSITY 1.//(KT*1000.*dt/(12.*PI*MONO_RAD*MONO_RAD*MONO_RAD))
#define FRICTION_COEFF 1.// 6pi*mono_rad*viscosity
#define MONO_MASS 1.
#define MONO_RAD 0.5//monomer diameter = 5nm (Lee and Liu, 2008)
#define MONO_RAD2 (MONO_RAD*MONO_RAD)
#define MONO_DIAM (2.*MONO_RAD)
#define MONO_DIAM2 (MONO_DIAM*MONO_DIAM)

#define OBSTACLE_RAD (10.*MONO_DIAM)
#define OBSTACLE_RAD2 (OBSTACLE_RAD*OBSTACLE_RAD)
#define OBSTACLE_MASS 1.
#define OBSTACLE_FRIC (20.*FRICTION_COEFF)


#define SIDE_BIND false
#define SPRINGCONST (cl_springconst)//((100.*KT)/(MONO_DIAM*MONO_DIAM))//soft repulsion coeff.
#define POLYMER_EXC_VOL_SPRING cl_polymer_exc_vol_spring
#define BOND_SPRING (cl_bond_spring)//((100.*KT)/(MONO_DIAM*MONO_DIAM))//bond potential coeff.
#define STIFFNESS (0.)//coefficient for bending potential
//coefficient for orientation potential
#define POLARIZABILITY (100.)//fixed

#define F_LOAD cl_f_load
#define LENGTH_CONTROLLED_LOAD cl_length_controlled_load
#define PIPETTE_STIFFNESS cl_pipette_stiffness
#define PIPETTE_VELOCITY cl_pipette_velocity

/*#define V0 0.
#define W0 0.*/

//parameters for geometric constraints for polymerization and interactions within filaments
//all but "orientation_cutoff" are taken from Lee and Liu (2008)
#define MIN_DIST_FOR_POLYM (0.9*MONO_DIAM)// = MONO_DIAM - 0.1*MONO_DIAM
#define MIN_DIST_FOR_POLYM_SQ (MIN_DIST_FOR_POLYM*MIN_DIST_FOR_POLYM)
#define FIL_ANGLE 0.
#define BRANCH_ANGLE (PI*70./180.)
//allows up to ~11.5degrees difference for filaments and ~68.8-71.2 for branches
#define INCOMING_ANGLE_CUTOFF 0.06//not actually an angle cutoff... a cut off for difference bet/ costheta and costheta0
#define ORIENTATION_CUTOFF 0.10//consider 0.1 but be aware that PE/force discontinuity is larger...
#define COS_F_ANGLE (cos(FIL_ANGLE))
#define COS_B_ANGLE (cos(BRANCH_ANGLE))


#define unif_rand() (RNUM1.get_double())//(ranf0()) // ((double) rand()/RAND_MAX) //uniform random number generator
#define ranf_std() (RNUM0.get_double())
#define gaussian_rand() (gaussian_inverf())//gaussian_std();

//ranges set to be 10*std dev.
#define sigma1sq (2.*TDIFF_COEFF*dt)
#define sigma2sq (2.*RDIFF_COEFF*dt)
#define sigma3sq (2.*TDIFF_COEFF/SHELL_RADIUS*dt)
#define sigma6sq (2.*RDIFF_COEFF*dt/DNA_DRAG_FACTOR)

#define RANGE_FACTOR1 (7.*sqrt(sigma1sq))//(10.*sqrt(sigma1sq))//gives max fluctuation
#define RANGE_FACTOR2 (10.*sqrt(sigma2sq))
#define RANGE_FACTOR3 (10.*sqrt(sigma3sq))
#define RANGE_FACTOR6 (10.*sqrt(sigma6sq))
                         /*Determines size of range of possible Gaussian random #'s generated.
                           Should be small enough so that we disallow crazy fluctuations.  Moreover,
                           making it small is advantageous because it speeds up the rejection method.
                           Making it too small with affect the accuracy of the Gaussian distribution, 
                           but 10*stddev is quite a large range...  P(10sigma) < 2e-22 * P(0)*/  



//#define RCUT2 (RCUT*RCUT)
#define invdt (1./dt)
#define invKT (1./KT)
#define TDIFF_COEFF (KT/(FRICTION_COEFF*MONO_MASS))
#define RDIFF_COEFF (0.75*TDIFF_COEFF/(MONO_RAD*MONO_RAD))
/*translational diffusion coefficient = kT/(6pi*eta*mono_rad)
rotational diffusion coefficient = kT/(8pi*eta*mono_rad^3)
See, for ex., Einstein's Investigations on the Theory of Brownian Motion*/
#define POLYMER_TDIFF_COEFF_FACTOR cl_polymer_tdiff_coeff_factor

//number of cells in each direction
#define NX (N[0])//25//14//9//((int)(LX/(4.1*MONO_DIAM)))
#define NY (N[1])//((int)(LY/(1.05*MONO_DIAM)))//25//9//((int)(LY/(4.1*MONO_DIAM)))
#define NZ (N[2])//((int)(LZ/(1.05*MONO_DIAM)))//25//9//((int)(LZ/(4.1*MONO_DIAM)))

//for neighbor_list/cell structure manipulations
#define LAST -1
#define EMPTY LAST

//for dynamic_monomer function calculate_1d_sep()
#define CURRENT false
#define PREVIOUS (!CURRENT)

//for dynamic_monomer function dynamical_interact()
#define FORM_BRANCH true
#define NOT_BRANCH (!FORM_BRANCH)

#define PI 3.141592653589793
#define TWOPI 6.283185307179586
#define HALFPI 1.5707963267949
#define ONE_OVER_HALF_SQRT_PI 0.28209479177387814
#define INV_PI 0.3183098861837907
#define SQRT_TWO 1.4142135623709505


/***********************************************************************************************************************************/
#define NUMBER_IN_POLYMER cl_number_of_parb_dna//0//360//0//100
#define NUMBER_OF_FILAMENTS 0
#define NUMBER_OF_SHELL_POLY_CONNECTIONS cl_number_of_shell_poly_connections
#define CROSSLINK_DENSITY cl_crosslink_density
#define NUMBER_OF_CROSSLINKS ((long)(0.5*CROSSLINK_DENSITY*NUMBER_IN_POLYMER))//((int)(6*NUMBER_OF_FILAMENTS)) //Number of crosslinks initially placed

#define POLYMER_CUT_PERCENTAGE cl_polymer_cut_percentage


//these only appear in crosslink_potential(), which is not used in current version of code 160928
#define CROSSLINK_SPRING (0.5*BOND_SPRING)
#define CROSSLINK_ATTRACTIVE_ONLY true
#define CROSSLINK_LENGTH (1.2*MONO_DIAM)
#define CROSSLINK_LENGTH2 (CROSSLINK_LENGTH*CROSSLINK_LENGTH)

#define PARB_STIFFNESS cl_parb_stiffness//100.

#define CROSSLINK true //crosslink the bundle

#define PARB_STICKS_NEAR_POLE false//if true, Parb sticks after it passes PARB_STICKING_PLANE
#define PARB_STICKING_PLANE (75.*MONO_DIAM)
#define PARB_DOES_NOT_BIND_AFTER_STICKING false



#define PRINT_CM_ONLY false
#define PRINT_RGYR true

//true -- side binding ring, false -- side binding sites
#define SIDE_BIND_RING false
//variables for side binding sites on ParA
#define RANDOM_SIDE_BIND_SITES false
#define NN_POLARIZATION_INT (0.)//strength of nearest neighbor polarization interaction
#define NEXT_NN_POLARIZATION_INT (0.*NN_POLARIZATION_INT)//strength of interaction between polarizations of next nearest neighbors
#define SIDE_ANGLE HALFPI//preferred angle between bond vector and polarization vector
#define COS_SIDE_ANGLE cos(SIDE_ANGLE)
#define COS_HALFPI_PLUS_B_ANGLE cos(BRANCH_ANGLE + HALFPI)//bind site will be oppposite branch
#define POL_ANGLE (PI*10./180.)//Alternatively, this could be defined as TWOPI/PITCH, if we define PITCH to be the number of monomers for the preferred polarization to rotate 2pi
#define COS_POL_ANGLE cos(POL_ANGLE)//preferred angle between nearest neighbor polarizations
#define COS_TWO_POL_ANGLE cos(2.*POL_ANGLE)//aka preferred angle between next nearest neighbor polarizations.

//variables for side binding sites on ParB
#define PARB_POLARIZABILITY 0.//(25.*KT)
#define PARB_NN_POLARIZATION_INT 0.//(50.*KT)




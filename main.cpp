#include "neighbor_list.h"
#include "random_numbers.cpp"
#include "RNG_taus.h"
#include "monomer_dynamics.cpp"
#include "update_system.cpp"
#include "binding.cpp"
#include "print_obstacle_pos.cpp"
#include "print_shell_spring_data.cpp"
#include "calculate_moduli_variables.cpp"
#include "create_shell.cpp"
#include "create_ellipsoid.cpp"
#include "create_cylinder.cpp"
#include "connect_shell_poly.cpp"
#include "quicksort.cpp"
#include "compress.cpp"
#include "filament_interactions.cpp"
#include "create_polymer.cpp"
#include "restart.cpp"
#include "initialize_sim.cpp"
#include "print_final_config.cpp"
#include "print_stresses.cpp"
#include "print_load_test.cpp"
#include "print_ylm.cpp"
#include "run_resize_box.cpp"
#include "write_state_to_pdb.cpp"

double shell_cm[DIMENSION];

//various lists of objects
vector<dynamic_monomer> mono_list;
dynamic_monomer central_mono;

vector<int> monolinklist;
vector<vector<vector<int> > > firstmonoincell;

vector<monomer_pair> pairs;
vector<monomer_pair> crosslinkpairs;


//some variables for rand num generators
int iseed0, iseed1, iseed2, iseed3, iseed4, iseed5, iseed6;
double gauss_prefact1, gauss_exp_const1, gauss_range1;
double gauss_prefact2, gauss_exp_const2, gauss_range2;
double gauss_prefact3, gauss_exp_const3, gauss_range3;
double gauss_prefact6, gauss_exp_const6, gauss_range6;
double gauss_prefact_std, gauss_exp_const_std, gauss_range_std;

//global variable saving the index of monomers that are indentation points if load force < 0
int indentation1, indentation2;

RNG_taus RNUM0, RNUM1;

double L[DIMENSION];
int N[DIMENSION];
double cell_length[DIMENSION], inv_cell_length[DIMENSION];

double center_mass_pos[DIMENSION];

int num_boundary;

vector<vector<bool> > attachment_matrix;
vector<bool> main_shell;
int tot_attachments;
vector<int> num_attachments;

int spring_mono_id;
double spring_att_point[DIMENSION];
double pipette_dist_from_center;

//rand # test
//long running_counter;

unsigned long long cl_numsteps;
int cl_trialnumber;
double cl_lx, cl_ly;
double cl_parb_stiffness;
double cl_f_load;
bool cl_length_controlled_load;
double cl_pipette_stiffness, cl_pipette_velocity;
double cl_bond_spring, cl_shell_bond_spring, cl_polymer_exc_vol_spring;
double cl_shell_radius;
int cl_rseed;
int cl_lower_bound_connectivity, cl_upper_bound_connectivity;
bool cl_thermal, cl_length_dependent_springs;
double cl_ldep_factor;
double cl_os_pressure;
bool cl_solid_interior;
double cl_springconst;
bool cl_springs_only, cl_extensional_springs_only;
bool cl_variable_bond_length, cl_modified_shell_exc_vol;
int cl_number_of_parb_dna, cl_num_shell_monos;
double cl_crosslink_spring_factor;
double cl_crosslink_density, cl_polymer_cut_percentage; 
int cl_number_of_shell_poly_connections;
int cl_num_load_monos;
double cl_load_frac;
double cl_polymer_mono_rad, cl_shell_mono_rad;
bool cl_print_lammps;
double cl_kt, cl_new_kt, cl_dt, cl_polymer_tdiff_coeff_factor;
bool cl_restart, cl_pdb;
bool cl_lowt_protocol, cl_box_resize, cl_extended_relaxation;
bool cl_cylinder;
double cl_cylinder_length;
bool cl_ellipsoid;
bool cl_calc_moduli;
bool cl_hysteresis;

bool cl_compress;
double TOP_WALLFORCE, BOTTOM_WALLFORCE, BOTTOM_WALLPOS, TOP_WALLPOS;

double saved_extension;

FILE *springfile;
char springfilename[128];
FILE *extfile;
char extfilename[128];
FILE *ylmfile;
char ylmname[128];
FILE *fourierfile;
char fouriername[128];
FILE *totforcefile;
char totforcename[128];
FILE *fluctuation_file;
char fluctuation_name[128];

void (*binding_func)(int, int);
void (*bending_func)(int);

int main(int argc, char **argv)
{
initialize_cl();

int option;
//remaining options are
//eJ
while( (option = getopt(argc, argv, "t:x:s:f:b:r:z:l:T:L:u:o:S:k:P:V:N:n:M:a:c:v:C:F:w:p:y:d:q:g:h:X:R:Z:A:B:W:Y:O:D:E:Q:m:G:H:i:K:j:U:I:")) != -1)
{
 switch(option)
 {
  case 't':
    cl_trialnumber = atoi(optarg); break;
  case 'n':
    cl_numsteps = atoi(optarg); break;
  case 'x':
    cl_lx = atof(optarg); break;
  case 's': 
    cl_parb_stiffness = atof(optarg); break;
  case 'f':
    cl_f_load = atof(optarg); break;
  case 'F':
    cl_load_frac = atof(optarg); break;
  case 'b':
    cl_bond_spring = atof(optarg); break;
  case 'a':
    cl_shell_bond_spring = atof(optarg); break;
  case 'c':
	cl_crosslink_density = atof(optarg); break;
  case 'C':
	cl_number_of_shell_poly_connections = atoi(optarg); break;
  case 'A':
	cl_polymer_cut_percentage = atof(optarg); break;
  case 'r':
    cl_shell_radius = atof(optarg);
    cl_ly = 4.05*cl_shell_radius;
//    if(cl_lx/NX < 1.0)
//	cl_lx = 1.05*NX; 
    if(4.0*cl_shell_radius > cl_lx)
	cl_lx = 4.05*cl_shell_radius;
    break;
  case 'o':
	cl_os_pressure = atof(optarg); break;
  case 'z':
    cl_rseed = atoi(optarg); break;
  case 'l':
        cl_lower_bound_connectivity = atoi(optarg); break;
  case 'u':
	cl_upper_bound_connectivity = atoi(optarg); break;
  case 'P':
	if(atoi(optarg) == 0)
	  cl_springs_only = false;
	break;
  case 'E':
	if(atoi(optarg) == 1)
	  cl_extensional_springs_only = true;
	break;
  case 'S':
	if(atoi(optarg) == 1)
	  cl_solid_interior = true;
	break;
  case 'k':
	cl_springconst = atof(optarg); break;
  case 'v':
	cl_polymer_exc_vol_spring = atof(optarg); break;
  case 'N':
	cl_number_of_parb_dna = atoi(optarg); break;
  case 'M':
	cl_num_shell_monos = atoi(optarg); break;
  case 'T':
	if(atoi(optarg) == 1)
	  cl_thermal = true;
	break;
  case 'L':
        if(atoi(optarg) == 1)
          cl_length_dependent_springs = true;
        break;
  case 'V':
	if(atoi(optarg) == 0)
	  cl_variable_bond_length = false;
	break;
  case 'w':
	cl_polymer_mono_rad = atof(optarg);
	break;
  case 'i':
	cl_shell_mono_rad = atof(optarg);
	break;
  case 'p':
	if(atoi(optarg) == 1)
	  cl_print_lammps = true;
	break;
  case 'y':
	cl_kt = atof(optarg);
	break;
  case 'd':
	cl_dt = atof(optarg);
	break;
  case 'q':
	if(atoi(optarg) == 1)
	  cl_length_controlled_load = true;
	break;
  case 'g':
	cl_pipette_stiffness = atof(optarg);
	break;
  case 'h':
	cl_pipette_velocity = atof(optarg);
	break;
  case 'X':
	if(atoi(optarg) == 1)
	  cl_modified_shell_exc_vol = true;
	break;
  case 'R':
	if(atoi(optarg) == 1)
	   cl_restart = true;
	break;
  case 'B':
	if(atoi(optarg) == 1)
	  cl_pdb = true;
	break;
  case 'Z':
	cl_polymer_tdiff_coeff_factor = atof(optarg);
	break;
  case 'W':
	if(atoi(optarg) == 0)
	  cl_lowt_protocol = false;
	break;
  case 'Y':
	if(atoi(optarg) == 0)
	  cl_box_resize = false;
	break;
  case 'O':
	if(atoi(optarg) == 1)
	  cl_cylinder = true;
	break;
  case 'D':
	cl_cylinder_length = atof(optarg); break;
  case 'Q':
	if(atoi(optarg) == 1)
	  cl_extended_relaxation = true;
	break;
  case 'm':
	cl_ldep_factor = atof(optarg); break;
  case 'G':
	if(atoi(optarg) == 1)
	 cl_ellipsoid = true;
	break;
  case 'H':
	if(atoi(optarg) == 1)
	  cl_calc_moduli = true;
	break;
  case 'K':
	if(atoi(optarg) == 1)
	  cl_hysteresis = true;
	break;
  case 'j':
	cl_new_kt = atof(optarg);
	break;
  case 'U':
	if(atoi(optarg) == 1) cl_compress = true;
	break;
  case 'I':
        cl_crosslink_spring_factor=atof(optarg);
        break;
  case '?':
	fprintf(stderr, "menu:\n");
	fprintf(stderr, "trial -t\n");
	fprintf(stderr, "dt -d\n");
	fprintf(stderr, "force -f\n");
	fprintf(stderr, "radius -r\n");
	fprintf(stderr, "stiffness -s\n");
	fprintf(stderr, "lx -x \n");
	fprintf(stderr, "seed -z\n");
        fprintf(stderr, "lower bound connectivity -l\n");
	fprintf(stderr, "num in polymer -N\n");
	fprintf(stderr, "num in shell -M\n");
	fprintf(stderr, "length dep spring consts -L\n");
	fprintf(stderr, "springs initialized at varied lengths -V\n");
	fprintf(stderr, "thermal noise -T\n");
	fprintf(stderr, "kT -y\n");
	fprintf(stderr, "polymer diff coeff factor -Z\n");
	fprintf(stderr, "restart -R\n");
	fprintf(stderr, "use modified excluded volume scheme -X\n");
	fprintf(stderr, "pipette stiffness -g\n");
	fprintf(stderr, "pipette velocity -h\n");
	fprintf(stderr, "and more...\n");
	fprintf(stderr, "Exiting\n");
	exit(1);
 }
}

fprintf(stderr, "Options chosen: t %i z %i r %g b %g f %g x %g\n", cl_trialnumber, cl_rseed, cl_shell_radius, cl_bond_spring, cl_f_load, cl_lx); 

//running_counter = 0;

   fprintf(stderr, "trial number %i will run for %llu steps.\n", TRIALNUMBER, NUMSTEPS);
   
   unsigned long long ii;
   long nn;

 
   if(RESTART)
     read_restart();
   else
     initialize_sim();  //initialize random number generator, sets id #'s for monomers, lots of other stuff 

int dummy;

//    char command[64];
//    sprintf(command, "mkdir output/data%6.6i", TRIALNUMBER);
//    system(command);

  char cpcommand[128];
  char mvcommand[128];
  sprintf(cpcommand, "cp /tmp/[frsy][oapl]*%6.6i output/", TRIALNUMBER);
  sprintf(mvcommand, "mv /tmp/[frsy][oapl]*%6.6i output/", TRIALNUMBER);

if((!RESTART) && LOWT_PROTOCOL)
{
bool LOW_TEMPERATURE_PROTOCOL;
if(KT < 0.5)
  LOW_TEMPERATURE_PROTOCOL = true;
else
  LOW_TEMPERATURE_PROTOCOL = false;

//movie
//LOW_TEMPERATURE_PROTOCOL = false;

if(LOW_TEMPERATURE_PROTOCOL || CALC_MODULI)//want transient for CALC_MODULI to resize box, set x0
{
vector<vector<double> > avgpos;
vector<double> temp_xvec(3,0.);
for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
{
 avgpos.push_back(temp_xvec);
}

	if(!SPRINGS_ONLY)
	{
	  for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
	  {
		mono_list[ii].set_tdiffusion_coeff(1.0);//decrease drag
		mono_list[ii].set_tdiff_stdev(sqrt(2.*TDIFF_COEFF*dt));//but keep the low temperature noise amplitude
	  }

	  unsigned long long num_transient_steps = 2000000;
          for(ii = 0; ii < num_transient_steps; ii++)
	  {
	   if(ii%(50*NUMSKIP) == 0)
	    fprintf(stderr, "low temp step %llu\n", ii);   
	
	   monomer_dynamics(false, ii);	
           filament_interactions();
           update_system(ii);

	   if(ii%(25*NUMSKIP) == 0)
	   {
		for(int ii2 = NUMBER_IN_POLYMER; ii2 < NUMBER_OF_MONOMERS; ii2++)
		{
		  for(int kk = 0; kk < DIMENSION; kk++)
		    avgpos[ii2][kk] += mono_list[ii2].get_prev_pos(kk) - shell_cm[kk];
		}

	   if(ii%(50*NUMSKIP) == 0)
	   {
	    run_resize_box();
	   }//if ii%(50*NUMSKIP
	   }//ii%25nsk
	  }//for(ii

          for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
          {
                mono_list[ii].set_tdiffusion_coeff(TDIFF_COEFF);//everything set to low temp now
		if(CALC_MODULI)
		if(ii >= NUMBER_IN_POLYMER)
		{
		  for(int kk = 0; kk < DIMENSION; kk++)
		  {
		   avgpos[ii][kk] *=  25.*(double)NUMSKIP / (double) num_transient_steps;
		   mono_list[ii].set_xbar(kk, avgpos[ii][kk]);
		  }
		}
	  }
        }//if(!SPRINGS_ONLY)
	else
	{
		fprintf(stderr, "low temp protocol not written for springs only system. continuing\n");
	}
}//if(LOW_TEMPERATURE_PROTOCOL || CALC_MODULI)
}//!RESTART && LOWT_PROTOCOL


if(HYSTERESIS)
{
   if(SPRINGS_ONLY)
   {  
     fprintf(stderr, "hysteresis currently only for normal runs. exiting\n");
     exit(1);
   }
   for(ii = 0; ii </*=*/ NUMSTEPS/2; ii++)
   {
      if(ii%(100*NUMSKIP) == 0)
      {
        fprintf(stderr, "hyst step number = %llu\n", ii);
      }

	monomer_dynamics(false, ii);
        filament_interactions();
        update_system(ii);

	if(BOX_RESIZE)
	if(ii % (NUMSKIP*100) == 0)
	run_resize_box();
/*
if(PDB)
if(ii%PDB_SKIP == 0)
    write_state_to_pdb(ii);
*/

     if((ii%500000 == 0) && (ii != 0))
      write_restart(ii);
   }//for loop over numstep/2 steps

cl_kt = cl_new_kt;
int n;
if(POLYMER_TDIFF_COEFF_FACTOR > 1.001){fprintf(stderr, "warning polymer tdiff factor > 1\n");}
for(n = 0; n < NUMBER_IN_POLYMER; n++)
{
  mono_list[n].set_tdiffusion_coeff(POLYMER_TDIFF_COEFF_FACTOR*TDIFF_COEFF);
}
for(n = NUMBER_IN_POLYMER; n < NUMBER_OF_MONOMERS; n++)
{
  mono_list[n].set_tdiffusion_coeff(TDIFF_COEFF);
}
}//if hysteresis

 
 
   for(ii = 0; ii </*=*/ NUMSTEPS; ii++)
   {
      if(ii%(100*NUMSKIP) == 0)
      {
        fprintf(stderr, "step number = %llu\n", ii); 
      }

if(!SPRINGS_ONLY)
{
	monomer_dynamics(false, ii);
}
else
{
	for(int qq = NUMBER_IN_POLYMER; qq < NUMBER_OF_MONOMERS; qq++)
	{
         if(mono_list[qq].get_load_mono())
         {
              (mono_list[qq]).load();
	 }
	
	 if(THERMAL)
	 {
              (mono_list[qq]).translational_brownian_move(ii);
	      if(SOLID_INTERIOR)
		 central_mono.translational_brownian_move(ii);	
         }

	 if(SOLID_INTERIOR)
	 {
		if(central_mono.calculate_distance_sq(&mono_list[qq], PREVIOUS) < SHELL_RADIUS*SHELL_RADIUS)
		{
		//	fprintf(stderr, "r_central= %g %g %g, r_mono= %g %g %g\n", central_mono.get_prev_pos(0), central_mono.get_prev_pos(1), central_mono.get_prev_pos(2), mono_list[qq].get_prev_pos(0), mono_list[qq].get_prev_pos(1), mono_list[qq].get_prev_pos(2));
			central_mono.sys_force(&mono_list[qq]);
		//	int wait;
		//	scanf("%i", &wait);
		}
	 }
	
	}
}//else, springs only

//fprintf(stderr, "r_central= %g %g %g\n", central_mono.get_prev_pos(0), central_mono.get_prev_pos(1), central_mono.get_prev_pos(2));
        filament_interactions();

	update_system(ii);


if(PDB)
/*
if(THERMAL && !SPRINGS_ONLY)
*/if(ii%PDB_SKIP == 0)
/*if(TRIALNUMBER % 22 == 0)
if((TRIALNUMBER%11 == 0) || (TRIALNUMBER > 999990) || ((TRIALNUMBER > 28435) && (TRIALNUMBER < 28442)) )
*/{
    write_state_to_pdb(ii);
}
/*#endif
*/
        if(ii%(NUMSKIP*2) == 0)
	{
          print_obstacle_pos(ii);
	  
	  if(CALC_MODULI)
		calculate_moduli_variables(ii);

	  #if PRINT_SHELL_SPRING
	  print_shell_spring_data(ii);
          if(ii > NUMSTEPS / 2)
          if(ii%(NUMSKIP*2000) == 0)
                print_stresses(ii);
	  #endif
	  #if (PRINT_FOURIER || PRINT_YLM)
	  print_ylm(ii);
	  #endif
	  #if LOAD_TEST
	  print_load_test(ii);
	  #endif
	}     


   if((ii%200000 == 0) && (ii != 0))
   {
    write_restart(ii);
    dummy= system(cpcommand);
   }


   }//for(....) numsteps
   
//if I uncomment this, make sure to remove setting of x0 after the transient
   //print_final_config();

   #if WRITE
   write_restart(ii);
   #endif
   fclose(extfile);
   #if PRINT_SHELL_SPRING
   fclose(springfile);
   #endif
   #if PRINT_YLM
   fclose(ylmfile);
   #endif
   #if PRINT_FOURIER
   fclose(fourierfile);
   #endif
   if(LENGTH_CONTROLLED_LOAD)
   {
     fclose(totforcefile);
   }
   if(CALC_MODULI)
	fclose(fluctuation_file);


   dummy= system(mvcommand);


fprintf(stderr, "trial completed.\n");
//   fprintf(stdout, "Press enter to end.\n");
   
//   char hold_on;
//   scanf("%c", &hold_on);
   return 0;
}




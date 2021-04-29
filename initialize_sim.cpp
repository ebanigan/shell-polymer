/*This function should be called at the beginning of each simulation. It initializes 
the random number generator and creates the monomer array.*/

void initialize_sim()
{

   long n, ii, jj;
   int kk, mm;
   int dd;

initialize_system_params();
initialize_files();
initialize_neighbor_list();

cl_num_load_monos = (int)(cl_num_shell_monos*cl_load_frac);

//initialize ParB-DNA subunits
for(n = 0; n < NUMBER_IN_POLYMER; n++)
{
  dynamic_monomer tempdmono(n, POLYMER_TDIFF_COEFF_FACTOR*TDIFF_COEFF);
  mono_list.push_back(tempdmono);
  mono_list[n].set_exc_vol_spring_const(POLYMER_EXC_VOL_SPRING);
  mono_list[n].set_radius(cl_polymer_mono_rad);
}


//now initialize monomers in shell
   for(n = NUMBER_IN_POLYMER; n < NUMBER_OF_MONOMERS; n++)
   {
      dynamic_monomer tempdmono(n, TDIFF_COEFF);
      mono_list.push_back(tempdmono);
      mono_list[n].set_radius(cl_shell_mono_rad);
   } // create rest of monomer list


//always initialize central mono?
/*if(SOLID_INTERIOR)
{*/
 dynamic_monomer tempdmono(0, TDIFF_COEFF/SHELL_RADIUS);
 central_mono = tempdmono;
 central_mono.set_radius(SHELL_RADIUS-MONO_RAD);//subtract mono rad b/c force subroutines add mono rad to overlap dist 
 for(kk = 0; kk < DIMENSION; kk++)
  central_mono.set_pos(kk, L[kk]*0.5);
 central_mono.set_prev_pos();
//}

//in the shell code OBVIOUSLY this must come next, or your code will crash!

if(NUMBER_OF_MONOMERS > num_shell_monos)
   create_polymer();

//fprintf(stderr, "after pol %i\n", pairs.size());
//#endif

//fprintf(stderr, "before shell %i\n", pairs.size());

if(CYLINDER)
{
  create_cylinder();
}
else
{
if(num_shell_monos > 0)
{
if(ELLIPSOID)
  create_ellipsoid();
else
  create_shell();
}
}
//fprintf(stderr, "after shell %i\n", pairs.size());

  update_system(0);

if(cl_print_lammps)
{
	print_lammps_file();
	fprintf(stderr, "Printed a file for LAMMPS. Exiting.\n");
	exit(1);
}

for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
  mono_list[ii].set_x0();

fprintf(stderr, "Initialization (part 1) complete. Beginning initial dynamics cycle.\n");

if(LOAD_TEST)
for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
        fprintf(stderr, "%li %g %g %g\n", ii, mono_list[ii].get_prev_pos(0), mono_list[ii].get_prev_pos(1), mono_list[ii].get_prev_pos(2));

initial_dynamics();

//let shell relax before connecting interior to shell
if(NUMBER_OF_SHELL_POLY_CONNECTIONS > 0)
        connect_shell_poly();

if(LENGTH_CONTROLLED_LOAD)
   set_pipette_position();

if(!COMPRESS)
{
if(num_shell_monos > 1)
  choose_load_monos();
}
else
{
 BOTTOM_WALLPOS = 1.0;
 TOP_WALLPOS = L[2] - 1.0;
 BOTTOM_WALLFORCE = 0.;
 TOP_WALLFORCE = 0.;
}

if(CALC_MODULI)
{
double max_y = 0.;
int max_y_mono = 0;
for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
//should make sure num load monos is set to 1 otherwise this will lead to some gibberish results
  if(mono_list[ii].get_load_mono())
	mono_list[ii].set_update_bool(false);
  else if(mono_list[ii].get_pos(1) > max_y)
  {
    max_y = mono_list[ii].get_pos(1);
    max_y_mono = ii;
  }
}
mono_list[max_y_mono].set_update_bool(false);
}

update_system(0);

/*
for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
{
	if(mono_list[ii].get_load_mono())
	 fprintf(stderr, "%i %g %g %g %g\n", ii, mono_list[ii].get_load_rest_length(), mono_list[ii].get_prev_pos(0), mono_list[ii].get_prev_pos(1), mono_list[ii].get_prev_pos(2));
}

*/

fprintf(stderr, "Initialization complete.\n");

/*
for(n = 0; n < NUMBER_OF_MONOMERS; n++)
{
fprintf(stderr, "m %i pos %g %g %g prev pos %g %g %g\n", n, mono_list[n].get_pos(0), mono_list[n].get_pos(1), mono_list[n].get_pos(2), mono_list[n].get_prev_pos(0), mono_list[n].get_prev_pos(1), mono_list[n].get_prev_pos(2));
}
for(n = 0; n < pairs.size(); n++)
{
	fprintf(stderr, "mono pair %i is %i %i with sep %g %g %g and bl %g\n", n, (*(pairs[n].first)).get_id(), (*(pairs[n].second)).get_id(), pairs[n].get_relative_position(0), pairs[n].get_relative_position(1), pairs[n].get_relative_position(2), pairs[n].get_bond_length());
}*/


}


void initial_dynamics()
{
long ii;
int kk;
for(ii = 0; ii < NUM_INIT_STEPS; ii++)
{
if(ii % (NUMSKIP*10) == 0)
  fprintf(stderr, "initial dynamics step = %li\n", ii);

monomer_dynamics(true, 0); 

filament_interactions();

update_system(0);
}


//recenter entire structure 
recenter();

  update_system(ii);

}


void initialize_system_params()
{
   int jj;

//initialize ranf(), the uniform generator   
   iseed0 = rseed+217;//TRIALNUMBER+217;
   iseed1 = iseed0 + 2;
   iseed2 = iseed1 + 2;
   iseed3 = iseed2 + 2;
   iseed4 = iseed3 + 2;
   iseed5 = iseed4 + 2;
   iseed6 = iseed5 + 2;
   RNUM0.set(iseed0+3);
   RNUM1.set(iseed0+5);

spring_mono_id = -1;

//initialize the gaussian generators
   gauss_prefact1 = 1./(sqrt(sigma1sq*TWOPI)); //prefactor, and maximum in probability distribution
   gauss_exp_const1 = -1./(2.*sigma1sq);
   gauss_range1 = 2.*RANGE_FACTOR1;//Want range to extend to least likely event on both sides of 0.
                                             //10sigma in each direction
   
   gauss_prefact2 = 1./(sqrt(sigma2sq*TWOPI));
   gauss_range2 = 2.*RANGE_FACTOR2;
   gauss_exp_const2 = -1./(2.*sigma2sq);
   
   gauss_prefact3 = 1./(sqrt(sigma3sq*TWOPI));
   gauss_range3 = 2.*RANGE_FACTOR3;
   gauss_exp_const3 = -1./(2.*sigma3sq);
   
   gauss_prefact6 = 1./(sqrt(sigma6sq*TWOPI));
   gauss_range6 = 2.*RANGE_FACTOR6;
   gauss_exp_const6 = -1./(2.*sigma6sq);


   gauss_prefact_std = 1./(sqrt(sigma1sq*TWOPI));
   gauss_range_std = 2.*RANGE_FACTOR1;
   gauss_exp_const_std = -1./(2.*sigma1sq);

   
//   fprintf(stderr, "Std dev for translation = %g\n", sqrt(sigma1sq));
//   fprintf(stderr, "Std dev for rotation = %g\n", sqrt(sigma2sq));
   

   if(BOX_RESIZE)
   {
        cl_lx = 4.05*SHELL_RADIUS;
        cl_ly = cl_lx;
   }

//Initialize L[] and N[] arrays so some procedures can be carried out in arbitrary dimension   
initialize_box();

pipette_dist_from_center = 0.;

 saved_extension = 0.;

//remnant of a trial of non-specific binding.
//a binding function could be inserted here instead of this empty function
//fprintf(stderr, "initializing null binding\n");
binding_func = &binding_null;

   if(POLYMER_STIFFNESS < 1.e-5)
	bending_func = &bending_null;
   else
	bending_func = &bending;



 BOTTOM_WALLPOS = 1.0;
 TOP_WALLPOS = L[2] - 1.0;
 BOTTOM_WALLFORCE = 0.;
 TOP_WALLFORCE = 0.;

}



void initialize_cl()
{
	cl_trialnumber = 999999;
	cl_lx = 50.;
	cl_polymer_stiffness = 0.;
	cl_f_load = 0.;
	cl_length_controlled_load = false;
	cl_pipette_stiffness = 100.;
	cl_pipette_velocity = 0.1;

	cl_shell_radius = 5.0;
        cl_ly = 4.05*SHELL_RADIUS;

	cl_bond_spring = 100.0;
        cl_shell_bond_spring = 1.0*BOND_SPRING;

	cl_rseed = 419;
        cl_lower_bound_connectivity = 3;
	cl_upper_bound_connectivity = 7;
	cl_thermal = false;
	cl_length_dependent_springs = false;
	cl_ldep_factor = 10.;
	cl_os_pressure = 0.0;
	cl_solid_interior = false;
	cl_springconst = 1.0;
	cl_polymer_exc_vol_spring = 100.0;
	cl_springs_only = true;
	cl_extensional_springs_only = false;
	cl_variable_bond_length = true;
	cl_modified_shell_exc_vol = false;
	cl_number_in_polymer = 0;
	cl_num_shell_monos = 200;
	cl_numsteps = (unsigned long long)20000000+1;
	cl_crosslink_density = 0.;
        cl_crosslink_spring_factor=1.0;
	cl_number_of_shell_poly_connections = 0;
	cl_polymer_cut_percentage = 0.;

	cl_load_frac = 0.1;
	cl_num_load_monos = ((int)(cl_load_frac*num_shell_monos));

	cl_polymer_mono_rad = 0.5;
	cl_shell_mono_rad = 0.5;
	cl_print_lammps = false;
	cl_kt = 1.0;
	cl_new_kt = cl_kt;
	cl_dt = 5.e-4;
	cl_pdb = false;
	cl_polymer_tdiff_coeff_factor = 1.;
	cl_hysteresis = false;

	cl_lowt_protocol = true;
	cl_box_resize = true;
	cl_extended_relaxation = false;

	cl_cylinder = false;
	cl_cylinder_length = SHELL_RADIUS;

	cl_calc_moduli = false;

	cl_compress = false;
}


void initialize_files()
{
#if PRINT_SHELL_SPRING
      sprintf(springfilename, "/tmp/spring_data%6.6i", TRIALNUMBER);
      springfile = fopen(springfilename, "w");
#endif

      sprintf(extfilename, "/tmp/radius_of_gyrationsq%6.6i", TRIALNUMBER);
      extfile = fopen(extfilename, "w");	

#if PRINT_YLM
      sprintf(ylmname, "/tmp/ylm_pert_data%6.6i", TRIALNUMBER);
      ylmfile = fopen(ylmname, "w");
#endif

#if PRINT_FOURIER
      sprintf(fouriername, "/tmp/fourier_pert_data%6.6i", TRIALNUMBER);
      fourierfile = fopen(fouriername, "w");
#endif

 if(LENGTH_CONTROLLED_LOAD)
 {
	sprintf(totforcename, "/tmp/forces%6.6i", TRIALNUMBER);
	totforcefile = fopen(totforcename, "w");
 }

if(CALC_MODULI)
{ 
 sprintf(fluctuation_name, "/tmp/fluctuations%6.6i", TRIALNUMBER);
 fluctuation_file = fopen(fluctuation_name, "w");
}
}



void choose_load_monos()
{
int ii,kk;
vector<double> dist_from_zero(num_shell_monos, 0.);
vector<double> dist_from_lx(num_shell_monos, 0.);
double sum_sq0, sum_sqlx;

//cycle through shell monos to choose
for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
   sum_sq0 = mono_list[ii].get_prev_pos(0)*mono_list[ii].get_prev_pos(0);
   sum_sqlx = (L[0]-mono_list[ii].get_prev_pos(0))*(L[0]-mono_list[ii].get_prev_pos(0)); 
   for(kk = 1; kk < DIMENSION; kk++)
   {
	sum_sq0 = sum_sq0 + (0.5*L[kk]-mono_list[ii].get_prev_pos(kk))*(0.5*L[kk]-mono_list[ii].get_prev_pos(kk));
	sum_sqlx = sum_sqlx + (0.5*L[kk]-mono_list[ii].get_prev_pos(kk))*(0.5*L[kk]-mono_list[ii].get_prev_pos(kk));
   }//why bother. why not.
   dist_from_zero[ii-NUMBER_IN_POLYMER] = sqrt(sum_sq0);
   dist_from_lx[ii-NUMBER_IN_POLYMER] = sqrt(sum_sqlx);
}

dist_from_zero = quicksort(dist_from_zero);
dist_from_lx = quicksort(dist_from_lx);

for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
   sum_sq0 = mono_list[ii].get_prev_pos(0)*mono_list[ii].get_prev_pos(0);
   sum_sqlx = (L[0]-mono_list[ii].get_prev_pos(0))*(L[0]-mono_list[ii].get_prev_pos(0));


   for(kk = 1; kk < DIMENSION; kk++)
   {
        sum_sq0 = sum_sq0 + (0.5*L[kk]-mono_list[ii].get_prev_pos(kk))*(0.5*L[kk]-mono_list[ii].get_prev_pos(kk));
        sum_sqlx = sum_sqlx + (0.5*L[kk]-mono_list[ii].get_prev_pos(kk))*(0.5*L[kk]-mono_list[ii].get_prev_pos(kk));
   }

   if(sqrt(sum_sq0) < dist_from_zero[NUM_LOAD_MONOS])
   {
	mono_list[ii].choose_load_mono(-1., pipette_dist_from_center);
	fprintf(stderr, "left load %i ", ii);
	fprintf(stderr, "%g %g %g", mono_list[ii].get_prev_pos(0), mono_list[ii].get_prev_pos(1), mono_list[ii].get_prev_pos(2));
        fprintf(stderr, "\n");

	if(F_LOAD < 0.)
	  indentation1 = ii;

   }
   else if((sqrt(sum_sqlx) < dist_from_lx[NUM_LOAD_MONOS]) && !LOAD_TEST)
   {
	mono_list[ii].choose_load_mono(1., pipette_dist_from_center);
	fprintf(stderr, "right load %i ", ii);
        fprintf(stderr, "%g %g %g", mono_list[ii].get_prev_pos(0), mono_list[ii].get_prev_pos(1), mono_list[ii].get_prev_pos(2));
        fprintf(stderr, "\n");

        if(F_LOAD < 0.)
          indentation2 = ii;

   }
   else if((sqrt(sum_sqlx) < dist_from_lx[1]) && LOAD_TEST)
   {
	spring_mono_id = ii;
	fprintf(stderr, "attached mono %i ", ii);
        fprintf(stderr, "%g %g %g", mono_list[ii].get_prev_pos(0), mono_list[ii].get_prev_pos(1), mono_list[ii].get_prev_pos(2));
	fprintf(stderr, "\n");
	for(kk = 0; kk < DIMENSION; kk++)
	   spring_att_point[kk] = mono_list[spring_mono_id].get_prev_pos(kk);
	spring_att_point[0] = spring_att_point[0] - MONO_RAD;
   }
}

//choose closest 10% to x=0 and LX
/*

make sure to put load into monodynamics
*/
//exit(1);
}




void print_lammps_file()
{
int ii, kk;
int num_atom_types = 2;//num_shell_monos+1;
int num_bond_types = 2;//pairs.size() - NUMBER_IN_POLYMER + 1;

char lammpsname[128];
FILE *lammpsfile;

sprintf(lammpsname, "output/lammps%6.6i", TRIALNUMBER);
lammpsfile= fopen(lammpsname, "w");

fprintf(lammpsfile, "LAMMPS configuration for r=%g N_poly=%i N_shell=%i\n", SHELL_RADIUS, NUMBER_IN_POLYMER, num_shell_monos);
fprintf(lammpsfile, "\n");

fprintf(lammpsfile, "%i atoms\n", NUMBER_OF_MONOMERS);
fprintf(lammpsfile, "%i bonds\n", int(pairs.size()));
fprintf(lammpsfile, "0 angles\n");
fprintf(lammpsfile, "0 dihedrals\n");
fprintf(lammpsfile, "0 impropers\n");
fprintf(lammpsfile, "\n");

fprintf(lammpsfile, "%i atom types\n", num_atom_types);
fprintf(lammpsfile, "%i bond types\n", num_bond_types);
//fprintf(lammpsfile, "0 angle types\n");
fprintf(lammpsfile, "\n");

fprintf(lammpsfile, "Masses\n");
fprintf(lammpsfile, "\n");
for(ii = 0; ii < num_atom_types; ii++)
 fprintf(lammpsfile, "%i 1\n", ii+1); 
fprintf(lammpsfile, "\n");

fprintf(lammpsfile, "Bond Coeffs\n");
fprintf(lammpsfile, "\n");
fprintf(lammpsfile, "%i %g\n", 1, BOND_SPRING); 
for(ii = 1; ii < num_bond_types; ii++)
 fprintf(lammpsfile, "%i %g\n", ii+1, SHELL_BOND_SPRING); 
fprintf(lammpsfile, "\n");

fprintf(lammpsfile, "Atoms\n");
fprintf(lammpsfile, "\n");
for(ii = 0; ii < NUMBER_IN_POLYMER; ii++)
{
	fprintf(lammpsfile, "%i ", ii+1);//monomer id	
	fprintf(lammpsfile, "1 ");//monomer type is polymer monomer
	for(kk = 0; kk < DIMENSION; kk++)
	 fprintf(lammpsfile, "%g ", mono_list[ii].get_pos(kk));//position
	fprintf(lammpsfile, "\n");
}
for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
        fprintf(lammpsfile, "%i ", ii+1);//monomer id   
	fprintf(lammpsfile, "2 ");//shell mono
//        fprintf(lammpsfile, "%i ", ii+2-NUMBER_IN_POLYMER);//monomer type is polymer monomer
	for(kk = 0; kk < DIMENSION; kk++)
         fprintf(lammpsfile, "%g ", mono_list[ii].get_pos(kk));//position
        fprintf(lammpsfile, "\n");
}
fprintf(lammpsfile, "\n");

fprintf(lammpsfile, "Bonds\n");
fprintf(lammpsfile, "\n");
for(ii = 0; ii < NUMBER_IN_POLYMER - 1; ii++)
{
	fprintf(lammpsfile, "%i ", ii+1);//bond id
	fprintf(lammpsfile, "%i ", 1);// bond type is polymer mono bond
	fprintf(lammpsfile, "%li %li ", (*(pairs[ii].first)).get_id()+1, (*(pairs[ii].second)).get_id()+1);//monomers involved in the bond
	fprintf(lammpsfile, "\n");
}
for(ii = NUMBER_IN_POLYMER - 1; ii < pairs.size(); ii++)
{
	fprintf(lammpsfile, "%i ", ii+1);//bond id
	fprintf(lammpsfile, "2 ");//bond type is shell mono bond..
//	fprintf(lammpsfile, "%i ", ii+2-NUMBER_IN_POLYMER);//bond type
	fprintf(lammpsfile, "%li %li ", (*(pairs[ii].first)).get_id()+1, (*(pairs[ii].second)).get_id()+1);//monomers involved in the bond
        fprintf(lammpsfile, "\n");
}
fprintf(lammpsfile, "\n");

fflush(lammpsfile);
fclose(lammpsfile);


}



void set_pipette_position()
{
int ii;
double min_shell = 9.*L[0];
double max_shell = -9.*L[0];
double monox;
double misalignment;

for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
        monox = mono_list[ii].get_prev_pos(0);
        if(monox < min_shell)
                min_shell = monox;
        else if(monox > max_shell)
                max_shell = monox;
}

misalignment = max_shell + min_shell - L[0];

for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
{
	mono_list[ii].move_1d(0, -0.5*misalignment);
	mono_list[ii].set_prev_pos();
}


pipette_dist_from_center = 0.5*(max_shell - min_shell);


/*
fprintf(stderr, "max %g min %g misalign %g pipette %g box %g\n", max_shell, min_shell, misalignment, pipette_dist_from_center, L[0]);

//find max and min again..
min_shell = 9.*L[0];
max_shell = -9.*L[0];
for(ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
        monox = mono_list[ii].get_prev_pos(0);
        if(monox < min_shell)
                min_shell = monox;
        else if(monox > max_shell)
                max_shell = monox;
}

fprintf(stderr, "new max %g new min %g pipette %g\n",  max_shell, min_shell, pipette_dist_from_center);
*/
}


void initialize_box()
{
   int jj;
   double interaction_dist;

   if(cl_shell_mono_rad > cl_polymer_mono_rad)
   {
        interaction_dist = 2.*cl_shell_mono_rad;
   }
   else
   {
        interaction_dist = 2.*cl_polymer_mono_rad;
   }


   double int_factor = 0.52;//gets multiplied by 2. in do loop, first value is thus 1.04
   unsigned long long numcells = 0;
   unsigned long long num_monos_2 = NUMBER_OF_MONOMERS*NUMBER_OF_MONOMERS;
   do
   {
    int_factor = 2.*int_factor;

    for(jj = 0; jj < DIMENSION; jj++)
    {
      if(jj == 0)
      {
        L[jj] = LX;
        N[jj] = (int)(LX / (int_factor*interaction_dist));
      }
      else if(jj == 1)
      {
        L[jj] = LY;
        N[jj] = (int)(LY / (int_factor*interaction_dist));
      }
      else if(jj == 2)
      {
        L[jj] = LZ;
        N[jj] = (int)(LZ / (int_factor*interaction_dist));
      }
      else
      {
        fprintf(stderr, "Dimension greater than 3. Exiting.\n");//4 dimensional cell crawling? no way!
        exit(1);
      }
      if(N[jj] < 3) // note that even if cells are smaller than they should be, monomer_dynamics cycles through all 3 directions properly since all cells are checked for interactions, even in the "finite box" case
      {
        fprintf(stderr, "correcting N[%i]<3\n", jj);
        N[jj] = 3;
      }
      cell_length[jj] = L[jj] / N[jj];
      inv_cell_length[jj] = 1./cell_length[jj];
    }//for loop


   numcells = N[0]*N[1]*N[2];
   }while((numcells > num_monos_2)&&(NUMBER_OF_MONOMERS > 27));


}


void recenter()
{
int ii, kk;

vector<double> low_shift(3,0.);
vector<double> high_shift(3,0.);
  
  for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
  {
        for(kk = 0; kk < DIMENSION; kk++)
        {
                mono_list[ii].set_pos(kk, mono_list[ii].get_pos(kk) - (shell_cm[kk] - 0.5*L[kk]));

		//recentering process is NOT guaranteed to keep all monos inside the box
		if(mono_list[ii].get_pos(kk) < 0.)
		{
		    if(mono_list[ii].get_pos(kk) < low_shift[kk])
		     low_shift[kk] = mono_list[ii].get_pos(kk);
		}
		else if(mono_list[ii].get_pos(kk) > L[kk])
		{
		    if(mono_list[ii].get_pos(kk) - L[kk] > high_shift[kk])
		     high_shift[kk] = mono_list[ii].get_pos(kk) - L[kk];
		}
        }
  }
  for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
  {
        for(kk = 0; kk < DIMENSION; kk++)
	  mono_list[ii].set_pos(kk, mono_list[ii].get_pos(kk) - 1.5*low_shift[kk] - 1.5*high_shift[kk]);

        mono_list[ii].set_prev_pos();
  }

  if(COMPRESS)
  {
	TOP_WALLPOS += shell_cm[kk] - 0.5*L[2];
	BOTTOM_WALLPOS += shell_cm[kk] - 0.5*L[2];
  }
}



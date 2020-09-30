#include "params.h"

extern int iseed0, iseed1, iseed2, iseed3, iseed4, iseed5, iseed6; //random # seed
extern double gauss_prefact1, gauss_exp_const1, gauss_range1;
extern double gauss_prefact2, gauss_exp_const2, gauss_range2;
extern double gauss_prefact3, gauss_exp_const3, gauss_range3;
extern double gauss_prefact6, gauss_exp_const6, gauss_range6;
extern double gauss_prefact_std, gauss_exp_const_std, gauss_range_std;
extern double L[DIMENSION];
extern int N[DIMENSION];
extern double cell_length[DIMENSION], inv_cell_length[DIMENSION];
extern double center_mass_pos[DIMENSION];

extern int indentation1, indentation2;

extern double LJ_MIN;
extern double LJ_MIN2;
extern double LJ_MIN6;
extern double LJ_MIN12;
extern double LJ_MIN14;
extern double LJ_MIN16;
extern double TETHER_RANGE;

//rand # test
//extern long running_counter;

extern int cl_trialnumber;
extern double cl_epsilon;
extern double cl_lx, cl_ly;
extern double cl_parb_stiffness;
extern double cl_f_load;
extern bool cl_length_controlled_load;
extern double cl_pipette_stiffness, cl_pipette_velocity;
extern double cl_bond_spring, cl_crosslink_spring_factor;
extern double cl_shell_radius;
extern int cl_rseed;
extern int cl_lower_bound_connectivity, cl_upper_bound_connectivity;
extern bool cl_springs_only, cl_thermal, cl_length_dependent_springs, cl_extensional_springs_only;
extern double cl_ldep_factor;
extern double cl_os_pressure;
extern bool cl_solid_interior;
extern double cl_springconst, cl_shell_bond_spring, cl_polymer_exc_vol_spring;
extern bool cl_variable_bond_length, cl_modified_shell_exc_vol;
extern int cl_number_of_parb_dna, cl_num_shell_monos;
extern unsigned long long cl_numsteps;
extern double cl_crosslink_density, cl_polymer_cut_percentage;
extern int cl_number_of_shell_poly_connections;
extern double shell_cm[DIMENSION];
extern int cl_num_load_monos;
extern double cl_load_frac;
extern double cl_polymer_mono_rad, cl_shell_mono_rad;
extern bool cl_print_lammps;
extern double cl_kt, cl_new_kt, cl_dt, cl_polymer_tdiff_coeff_factor;
extern bool cl_restart, cl_pdb;

extern bool cl_lowt_protocol, cl_box_resize, cl_extended_relaxation, cl_calc_moduli;

extern bool cl_cylinder;
extern double cl_cylinder_length;
extern bool cl_ellipsoid;
extern bool cl_hysteresis;

extern int spring_mono_id;
extern double spring_att_point[DIMENSION];

extern double pipette_dist_from_center;

extern int num_boundary;

extern bool cl_compress;
extern double TOP_WALLFORCE, BOTTOM_WALLFORCE, BOTTOM_WALLPOS, TOP_WALLPOS;

using namespace std;
extern vector<vector<bool> > attachment_matrix;
extern vector<bool> main_shell;
extern int tot_attachments;
extern vector<int> num_attachments;

extern double saved_extension;

extern FILE *springfile;
extern char springfilename[128];
extern FILE *extfile;
extern char extfilename[128];
extern FILE *ylmfile;
extern char ylmname[128];
extern FILE *fourierfile;
extern char fouriername[128];
extern FILE *totforcefile;
extern char totforcename[128];
extern FILE *fluctuation_file;
extern char fluctuation_name[128];

extern void (*binding_func)(int, int);
extern void (*bending_func)(int);


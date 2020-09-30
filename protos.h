#include<vector>
using namespace std;

extern double ranf0(); //rand # gen (uniform, 0 to 1)
extern double ranf1();
extern double ranf2();
extern double ranf3();
extern double ranf4();
extern double ranf5();
extern double ranf6();
extern double gaussian1();//gaussian rand # gen. for translational brownian motion of monomers
extern double gaussian2();//rot brownian for monos
extern double gaussian3();
extern double gaussian6();
extern double gaussian_std();
extern double gaussian_inverf(unsigned long long);


extern void write_restart(unsigned long long);
extern void read_restart();

extern void initialize_sim(); //initializes rand # gen, system, etc.
extern void initial_dynamics();
extern void initialize_system_params();
extern void initialize_files();
extern void initialize_box();
extern void recenter();
extern void run_resize_box();
extern void resize_box(double, double, double);
extern void visit(int);
extern void add_bond(int, int);
extern void set_pipette_position();

extern void update_system(long);

//extern void actin_event_loop2(); //runs through possible actin events
extern void filament_interactions();

void write_state_to_pdb(unsigned long long);


extern void create_polymer();
extern void create_cylinder();
extern void connect_shell_poly();
extern void relax_crosslinks(double);
extern void update_att_times();

extern void actin_crosslink();
extern void remove_extra_actin();


extern void binding(int, int);
extern void binding_null(int, int);
extern void print_stresses(unsigned long long);
extern void print_obstacle_pos(unsigned long long);
extern void print_shell_spring_data(unsigned long long);
extern void print_load_test(long);
extern void print_final_config();
extern void monomer_dynamics(bool, unsigned long long);
extern void compress(int);
extern void initialize_cl();
extern void print_lammps_file();
extern void calculate_moduli_variables(unsigned long long);

extern void create_shell();
extern double choose_shell_spring_const(double);
extern void choose_shell_attachments();
extern void choose_load_monos();
vector<double> quicksort(vector<double>);

extern void create_ellipsoid();

int get_firstmonoincell(int, int, int);
int get_monolinklist(int);
void initialize_neighbor_list();
void update_mono_list();


extern void print_ylm(unsigned long long);
extern double y00(double, double, double, double);
extern double y10(double, double, double, double);
extern double y20(double, double, double, double);
extern double y30(double, double, double, double);
extern double y40(double, double, double, double);

extern double re_y11(double, double, double, double);
extern double im_y11(double, double, double, double);
extern double re_y21(double, double, double, double);
extern double im_y21(double, double, double, double);
extern double re_y22(double, double, double, double);
extern double im_y22(double, double, double, double);

extern double re_exp_nphi(double, double, double, int);
extern double im_exp_nphi(double, double, double, int);


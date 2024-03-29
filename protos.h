#include<vector>
//function prototypes

using namespace std;

extern double ranf0(); //rand # gen (uniform, 0 to 1)
extern double ranf2();
extern double ranf5();
extern double ranf6();
extern double gaussian2();//rot brownian for monos
extern double gaussian6();
extern double gaussian_std();
extern double gaussian_inverf(unsigned long long);

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
extern void polymer_interactions();

void write_state_to_pdb(unsigned long long);

extern void create_polymer();
extern void create_cylinder();
extern void connect_shell_poly();
extern void relax_crosslinks(double);

extern void binding(int, int);
extern void binding_null(int, int);
extern void print_stresses(unsigned long long);
extern void print_ext(unsigned long long);
extern void print_shell_spring_data(unsigned long long);
extern void print_load_test(long);
extern void monomer_dynamics(bool, unsigned long long);
extern void compress(int);
extern void initialize_cl();

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

extern void print_fourier(unsigned long long);
extern double re_exp_nphi(double, double, double, int);
extern double im_exp_nphi(double, double, double, int);


#define _CRT_SECURE_NO_WARNINGS
#ifndef _GLOBAL_H_
#define _GLOBAL_H_ 1

#define INF 1.0e14
#define EPS 1.0e-10
#define E 2.71828182845905
#define PI 3.14159265358979
#define GNUPLOT_COMMAND "gnuplot -persist"

typedef struct 
{
    int rank;
    double constr_violation;
    double *xreal;
    int **gene;
    double *xbin;
    double *obj;
    double *constr;
    double realvar_crowd_dist;
    double binvar_crowd_dist;
    double obj_crowd_dist;
    double crowd_dist;
} individual;

typedef struct
{
    individual *ind;
} population;

typedef struct lists
{
    int index;
    struct lists *parent;
    struct lists *child;
} list;

extern int nreal;
extern int nbin;
extern int nobj;
extern int ncon;
extern int popsize;
extern double pcross_real;
extern double pcross_bin;
extern double pmut_real;
extern double pmut_bin;
extern double eta_c;
extern double eta_m;
extern int ngen;
extern int nbinmut;
extern int nrealmut;
extern int nbincross;
extern int nrealcross;
extern int *nbits;
extern double *min_realvar;
extern double *max_realvar;
extern double *min_binvar;
extern double *max_binvar;
extern double *min_obj;
extern double *max_obj;
extern double *epsilon;
extern int bitlength;
extern double delta;
extern int mate;
extern int input_type;
extern int run_mode;
extern int var_option;
extern int obj_option;
extern int frequency;
extern int var1;
extern int obj1;
extern int obj2;
extern int obj3;
extern int angle1;
extern int angle2;
extern int choice;

void allocate_memory_pop (population *, int);
void allocate_memory_ind (individual *);
void deallocate_memory_pop (population *, int);
void deallocate_memory_ind (individual *);

void crossover (individual *, individual *, individual *, individual *);
void realcross (individual *, individual *, individual *, individual *);
void bincross (individual *, individual *, individual *, individual *);

void assign_crowding_distance_list (population *, list *, int);
void assign_crowding_distance_list_obj (population *, list *, int);
void assign_crowding_distance_list_realvar (population *, list *, int);
void assign_crowding_distance_list_binvar (population *, list *, int);
void assign_crowding_distance_indices (population *, int, int);
void assign_crowding_distance_indices_obj (population *, int, int);
void assign_crowding_distance_indices_realvar (population *, int, int);
void assign_crowding_distance_indices_binvar (population *, int, int);

void decode_pop(population *);
void decode_ind(individual *);

void onthefly_display (population *, FILE *, int);

int check_loose_dominance (individual *, individual *);
int check_dominance (individual *, individual *);
void define_epsilon (population *, int , double *);

void evaluate_pop (population *);
void evaluate_ind (individual *);

void fill_nondominated_sort (population *, population *);
void crowding_fill (population *, population *, int , int , list *);

void ga_optimizer(void);

void input(void);

void initialize_pop (population *);
void initialize_random_pop (population *);
void initialize_latin_pop (population *);
void initialize_file_pop (population *);

void insert (list *, int);
list* del (list *);

void merge(population *, population *, population *);
void copy_ind (individual *, individual *);

void mutation_pop (population *);
void mutation_ind (individual *);
void bin_mutate_ind (individual *);
void real_mutate_ind (individual *);

void output (FILE *, FILE *, FILE *, FILE *, FILE *);

void pdefinit(void);

void test_problem (double *, double *, int **, double *, double *);

void assign_rank_and_crowding_distance (population *);

void report_pop (population *, FILE *);
void report_feasible (population *, FILE *);
void report_var (population *, FILE *);
void report_obj (population *, FILE *);

void quicksort_front_obj(population *, int , int *, int);
void q_sort_front_obj(population *, int , int *, int , int );
void quicksort_front_realvar(population *, int , int *, int );
void q_sort_front_realvar(population *, int , int *, int , int );
void quicksort_front_binvar(population *, int , int *, int );
void q_sort_front_binvar(population *, int , int *, int , int );
void quicksort_dist(population *, int *, int );
void q_sort_dist(population *, int *, int , int );

void selection (population *, population *);
void restricted_selection (population *, population *);
individual* tournament (individual *, individual *);
list* search_nearest (population *, list *, list *, double *, double *, double *, double *);
list* search_second_nearest (population *, list *, list *, list *, double *, double *, double *, double *);
double calc_distance (population *, list *, list *, double *, double *, double *, double *);

double maximum (double , double);
double minimum (double , double);
double square (double);
double modulus (double *, int);
double dot (double *, double *, int);
double mean (double *, int);

#endif

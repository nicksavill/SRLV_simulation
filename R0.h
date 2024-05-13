#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>

#define min(A, B) ((A) < (B) ? (A) : (B))
#define max(A, B) ((A) > (B) ? (A) : (B))

enum sero_state {Negative, Positive};
enum state {Removed, Susceptible, Exposed, Infectious, Infected, Any, Resistant};
enum model_type{Next_Gen_Matrix, Deterministic, Stochastic, Parameters};
enum month {April, May, June, July, August, September, October, November, Decemeber, January, February, March};
enum inf_fn {Constant, Linear, Exponential, Four_yos};
enum strategy {Random, Oldest, Poor};

typedef struct Params {
    enum inf_fn infectious_func;
    enum month weaning_month;
    enum month sale_month;
    enum month mating_month;
    enum strategy sale_strategy;
    int num_ibm_sims;
    int max_age;
    int sim_years;
    int n_age_classes;
    int housing_period;
    int age_first_mating;
    int tsi_removal;
    double flock_size;
    double all_ewes_flock_size;
    double pens;
    double inter_pen_contact_weight;
    double beta_field;
    double beta_housed;
    double L_mean;
    double L_sd;
    double L_a;
    double L_b;
    double S_mean;
    double S_sd;
    double S_a;
    double S_b;
    double prob_dam_to_lamb;
    double prop_resistant;
    double gr_tsi_grad;
    double (*inf_fn)(int, struct Params*);
    double inf_fn_scale;
    double ***(*contact_weight_func)(struct Params*);
    char text[200];
    char file[200];
    gsl_rng *stream;
} params_t;

typedef struct {
    int sim_years, max_age;
    double **infected_by_cohort;
    double **positive_by_age;
    double **susceptible_by_age;
    double **infected_by_age;
    double **infectious_by_age;
    double **resistant_by_age;
    double **infectiousness_by_age;
    double **foi_by_age;
    double **total_by_age;
    double **tsi;
    double *positive;
    double *susceptible;
    double *infected;
    double *resistant;
    double *infectious;
    double *total;
    double *flock_size;
    double *positive_flock;
    double *infected_flock;
    double *infectious_flock;
    double *mean_infectious_age;
    double *seroconversions;
    double *prop_pos;
    int init_time;
    double R0, generation_time, doubling_time;
} result_t;

void set_seedRNG(int seed, params_t *p);
int uniform_intRNG(int a, int b, const gsl_rng *stream);
int binomialRNG(int n, double p, const gsl_rng *stream);
double gammaRNG(double a, double b, const gsl_rng *stream);
int* intvector(size_t size);
double* doublevector(size_t size);
double** doublematrix(size_t sizex, size_t sizey);
double*** doublematrix3(size_t sizex, size_t sizey, size_t sizez);
double size(int month, int start_age, int end_age, int n_age_classes, double *cohort);
int get_group(int age, int group_max_age[], int num_groups);
double* prob_infectious(params_t *params);
double* natural_deaths(int n_age_classes);
double* monthly_infection_rate(int housing_period, double beta_field, double beta_housed);
double*** contact_weights(params_t *p);
int* group_maximum_ages(params_t *params);
double* cohort_size(params_t *params, double *removal_rate);
double* cumulative_removal_rate(int n_age_classes, double *removal_rate);
double mean_infected_life_expectancy(params_t *params, double *uninfected_cohort_size, double *Lambda);
result_t* init_results(int sim_years, int max_age, int n_age_classes);
void result_free(result_t *r, enum model_type type);
void update_parameters(params_t *p);
void output_result(FILE *fp, result_t *r, params_t *p, enum model_type type);
double*** init_contact_weights(int max_age);
void run(params_t *p, params_t *init_p);
void run_single(params_t *p);
void run_multiple(params_t *p);
void get_retain_probabilities(double *alpha, double *uninfected_cohort_size, double *lambs_per_dam, params_t *p);
double** litter_dist_by_age_data(int age_first_mating, int max_age);
void print_contact_weights(double ***contact_weights, enum month m, params_t *p);

extern double constant_inf(int, params_t*);
extern double linear_inf(int, params_t*);
extern double exponential_inf(int, params_t*);
extern double four_yos_inf(int, params_t*);

result_t* Stochastic_Sim(params_t*, void**);
result_t* Deterministic_Sim(params_t*);
result_t* Next_Gen_Matrix_R0(params_t*);

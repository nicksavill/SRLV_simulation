#include "R0.h"

int* intvector(size_t size) {
    int *m = (int*) calloc(size, sizeof(int));
    if (!m) {
        fprintf(stderr, "Can't allocate %ld memory in intvector()", size);
        exit(0);
    }
    return m;
}

double* doublevector(size_t size) {
    double *m = (double*) calloc(size, sizeof(double));
    if (!m) {
        fprintf(stderr, "Can't allocate %ld memory in doublevector()", size);
        exit(0);
    }
    return m;
}

double** doublematrix(size_t sizex, size_t sizey) {
    size_t i;
    double **m = (double **) calloc(sizex, sizeof(double*));
    if (!m) {
        fprintf(stderr, "Can't allocate %ld memory in doublematrix()", sizex*sizey);
        exit(0);
    }
    m[0] = (double *) calloc(sizex * sizey, sizeof(double));
    if (!m[0]) {
        fprintf(stderr, "Can't allocate %ld memory in doublematrix()", sizex*sizey);
        exit(0);
    }
    for (i = 1; i < sizex; i++)
        m[i] = m[i-1] + sizey;
    return m;
}

double*** doublematrix3(size_t sizex, size_t sizey, size_t sizez) {
    size_t i, j;
    double ***m;
    double *lastaddress;

    m = (double ***) calloc(sizex, sizeof(double**));
    if (!m) {
        fprintf(stderr, "Can't allocate %ld memory in doublematrix()", sizex*sizey*sizez);
        exit(0);
    }
    m[0] = (double **) calloc(sizex * sizey, sizeof(double*));
    if (!m[0]) {
        fprintf(stderr, "Can't allocate %ld memory in doublematrix()", sizex*sizey*sizez);
        exit(0);
    }
    for (i = 1; i < sizex; i++)
        m[i] = m[i-1] + sizey;

    m[0][0] = (double *) calloc(sizex * sizey * sizez, sizeof(double));
    if (!m[0][0]) {
        fprintf(stderr, "Can't allocate %ld memory in doublematrix()", sizex*sizey*sizez);
        exit(0);
    }
    for (j = 1; j < sizey; j++)
        m[0][j] = m[0][j-1] + sizez;

    lastaddress = m[0][sizey-1];
    for (i = 1; i < sizex; i++)
        for (j = 0; j < sizey; j++) {
            m[i][j] = lastaddress + sizez;
            lastaddress = m[i][j];
        }
    return m;
}

double size(int month, int start_age, int end_age, int n_age_classes, double *uninfected_cohort_size) {
    // get the number of ewes from ages 12*start_age+month to 12*end_age+month
    int y, age;
    double sum = 0;
    for (y = start_age; y <= end_age; y++) {
        age = 12*y+month;
        if (age < n_age_classes)
            sum += uninfected_cohort_size[age];
    }
    return sum;
}

double* prob_infectious(params_t *p) {
    // Probability a ewe is infectious for each month after infection
    // Assumed to be Gamma distributed with mean L_mean and
    // standard deviation L_sd
    int a;
    double *p_infectious = doublevector(p->n_age_classes);

    for (a = 0; a < p->n_age_classes; a++) {
        // p_infectious[a] = 1;
        p_infectious[a] = gsl_cdf_gamma_P(a, p->L_a, p->L_b);
        // printf("%d %.30lf\n", a, p_infectious[a]);
    }
    // exit(0);
    return p_infectious;
}

double* natural_deaths(int n_age_classes) {
    int a;
    double p_mort, p_cull;
    double *removal_rate = doublevector(n_age_classes);

    // age-specific mortality and cull probabilities from Andrew's data
    double prob_mortality[] = {0.001, 0.012, 0.027, 0.041, 0.019, 0.021, 0.034};
    double prob_culled[] = {0, 0.046, 0.065, 0.093, 0.134, 0.196, 0.295};
    // assume no mortality or culling
    // double prob_mortality[] = {0, 0, 0, 0, 0, 0, 0};
    // double prob_culled[] = {0, 0, 0, 0, 0, 0, 0};

    // vector of removal rates
    for (a = 0; a < n_age_classes; a++) {
        p_mort = prob_mortality[a/12];
        p_cull = prob_culled[a/12];
        removal_rate[a] = -(log(1. - p_mort)+log(1. - p_cull))/12.;
    }
    return removal_rate;
}

double* monthly_infection_rate(int housing_period, double beta_housed, double beta_field) {
    // infection rate per infectious ewe each month
    int m;
    double *infection_rate = doublevector(12);

    if (housing_period > 12) {
        fprintf(stderr, "housing_period should be at most 12 months");
        exit(0);
    }

    for (m = April; m <= March - housing_period; m++)
        infection_rate[m] = beta_field;
    for (; m <= March; m++)
        infection_rate[m] = beta_housed;
    return infection_rate;
}

double* cohort_size(params_t *p, double *removal_rate) {
    // Calculate the number of individuals in an uninfected cohort through its lifespan
    // The mating month is the month the flock is counted and should equal flock_size
    int n_age_classes = p->n_age_classes;
    int max_age = p->max_age;
    int sale_month = p->sale_month;
    int mating_month = p->mating_month;
    int age_first_mating = p->age_first_mating;
    double flock_size = p->flock_size;

    int a, i;
    double sum;
    double *uninfected_cohort_size = doublevector(n_age_classes);

    // Set the sale month flock size to 1 (could be anything)
    uninfected_cohort_size[sale_month] = 1;

    // cohort size at age a after sale age until last month in flock
    for (a = sale_month+1; a < n_age_classes; a++)
        uninfected_cohort_size[a] = uninfected_cohort_size[a-1]*exp(-removal_rate[a-1]);

    // sum size of all existent cohorts at the mating month
    // the number of mated ewes in the mating month equals flock_size
    // the number of mated ewes does not include lambs and
    // ewes younger that 12*age_first_mating+mating_month (19 months)
    sum = size(mating_month, age_first_mating, max_age, n_age_classes, uninfected_cohort_size);

    // scale the cohort sizes with age 
    for (a = sale_month; a < n_age_classes; a++)
        uninfected_cohort_size[a] *= flock_size/sum;

    double **litter_dist_by_age = litter_dist_by_age_data(age_first_mating, max_age);
    uninfected_cohort_size[April] = 0;
    for (a = age_first_mating+1; a <= max_age; a++)
        for (i = 0; i < 5; i++)
            uninfected_cohort_size[April] += uninfected_cohort_size[12*a + April] * i * litter_dist_by_age[a][i];

    // dams = size(April, age_first_mating+1, max_age, n_age_classes, uninfected_cohort_size);
    // printf("dams at birth = %.1f\n", dams);

    // to get zero sold lambs, back-calculate number of lambs and set lambs_per_dam
    // for (a = sale_month-1; a >= 0; a--)
    //     uninfected_cohort_size[a] = uninfected_cohort_size[a+1]*exp(removal_rate[a+1]);
    // p->lambs_per_dam = uninfected_cohort_size[0]/dams;
    
    // number of lambs: cohort size at age 0
    // uninfected_cohort_size[April] = p->lambs_per_dam * dams;

    // cohort size at age a from birth to just before sale age
    for (a = May; a < sale_month; a++)
        uninfected_cohort_size[a] = uninfected_cohort_size[a-1]*exp(-removal_rate[a-1]);

    // for (a = 0; a < n_age_classes; a++)
    //     printf("%d %f\n", a, uninfected_cohort_size[a]);
    // exit(0);
    p->all_ewes_flock_size = size(April, age_first_mating, max_age, n_age_classes, uninfected_cohort_size);

    free(litter_dist_by_age[0]);
    free(litter_dist_by_age);
    
    return uninfected_cohort_size;
}

double* cumulative_removal_rate(int n_age_classes, double *removal_rate) {
    // cumulative rate of removal
    int a;
    double *Lambda = doublevector(n_age_classes);

    Lambda[0] = 0;
    for (a = 1; a < n_age_classes; a++)
        Lambda[a] = Lambda[a-1]+removal_rate[a-1];
    return Lambda;
}

double mean_infected_life_expectancy(params_t *p, double *uninfected_cohort_size, double *Lambda) {
    // calculate weigthed mean life expectancy post infection
    int a, t;
    int n_age_classes = p->n_age_classes;
    double life_expectancy[n_age_classes];
    double size_sum = 0, size_life_sum = 0;

    for (a = 0; a < n_age_classes; a++) {
        life_expectancy[a] = 0;
        // expected life expectancy from age of infection a
        // life expectancy since becoming infected is the sum of the 
        // probabilities of begin alive in age class from the age of becoming 
        // infected until removal
        for (t = a; t < n_age_classes; t++)
            life_expectancy[a] += exp(-Lambda[t]+Lambda[a]);
        // weight life expectancy by number of ewes in cohort at age a
        size_life_sum += uninfected_cohort_size[a] * life_expectancy[a];
        size_sum += uninfected_cohort_size[a];
    }
    return size_life_sum/size_sum;
}

result_t* init_results(int sim_years, int max_age, int n_age_classes) {
    result_t *result = (result_t*) malloc(sizeof(result_t));
    if (sim_years) {
        result->sim_years = sim_years;
        result->max_age = max_age;
        result->infected_by_cohort = doublematrix(sim_years+max_age+1, (max_age+1)*12);
        result->positive_by_age = doublematrix(12*sim_years, max_age+1);
        result->susceptible_by_age = doublematrix(12*sim_years, max_age+1);
        result->infected_by_age = doublematrix(12*sim_years, max_age+1);
        result->infectious_by_age = doublematrix(12*sim_years, max_age+1);
        result->resistant_by_age = doublematrix(12*sim_years, max_age+1);
        result->infectiousness_by_age = doublematrix(12*sim_years, max_age+1);
        result->foi_by_age = doublematrix(12*sim_years, max_age+1);
        result->total_by_age = doublematrix(12*sim_years, max_age+1);
        result->tsi = doublematrix(12*sim_years, max_age+1);
        result->positive = doublevector(12*sim_years);
        result->susceptible = doublevector(12*sim_years);
        result->infected = doublevector(12*sim_years);
        result->infectious = doublevector(12*sim_years);
        result->resistant = doublevector(12*sim_years);
        result->total = doublevector(12*sim_years);
        result->positive_flock = doublevector(12*sim_years);
        result->infected_flock = doublevector(12*sim_years);
        result->infectious_flock = doublevector(12*sim_years);
        result->flock_size = doublevector(12*sim_years);
        result->mean_infectious_age = doublevector(12*sim_years);
        result->seroconversions = doublevector(12*sim_years);
        result->prop_pos = doublevector(12*sim_years);
    }
    return result;
}

void result_free(result_t *r, enum model_type type) {
    if (type != Next_Gen_Matrix) {
        free(r->infected_by_cohort[0]);
        free(r->positive_by_age[0]);
        free(r->susceptible_by_age[0]);
        free(r->infected_by_age[0]);
        free(r->total_by_age[0]);
        free(r->infected_by_cohort);
        free(r->positive_by_age);
        free(r->susceptible_by_age);
        free(r->infected_by_age);
        free(r->infectious_by_age);
        free(r->resistant_by_age);
        free(r->infectiousness_by_age);
        free(r->foi_by_age);
        free(r->total_by_age);
        free(r->tsi);
        free(r->positive);
        free(r->susceptible);
        free(r->infected);
        free(r->infectious);
        free(r->resistant);
        free(r->total);
        free(r->positive_flock);
        free(r->infected_flock);
        free(r->infectious_flock);
        free(r->flock_size);
        free(r->mean_infectious_age);
        free(r->seroconversions);
        free(r->prop_pos);
    }
    free(r);
}

void update_parameters(params_t *p) {
    if (p->contact_weight_func == NULL) {
        fprintf(stderr, "set contact weight function\n");
        exit(0);
    }

    // total number of age classes (maximum number of months a ewe can stay in the flock)
    p->n_age_classes = 12*p->max_age + p->sale_month;

    // gamma distribution parameters for latent period
    p->L_a = pow(p->L_mean/p->L_sd, 2.);
    p->L_b = p->L_sd*p->L_sd/p->L_mean;

    // gamma distribution parameters for seroconversion period
    p->S_a = pow(p->S_mean/p->S_sd, 2.);
    p->S_b = p->S_sd*p->S_sd/p->S_mean;

    if (p->infectious_func == Constant)
        p->inf_fn = &constant_inf;
    else if (p->infectious_func == Linear)
        p->inf_fn = &linear_inf;
    else if (p->infectious_func == Exponential)
        p->inf_fn = &exponential_inf;
    else if (p->infectious_func == Four_yos)
        p->inf_fn = &four_yos_inf;
    else {
        fprintf(stderr, "infectiousness function undefined\n");
        exit(0);
    }

  // scale of infectious function
    int a;

    p->inf_fn_scale = 0;
    for (a = 0; a < p->n_age_classes; a++)
        p->inf_fn_scale += (*p->inf_fn)(a, p);
    p->inf_fn_scale = p->n_age_classes / p->inf_fn_scale;
    // printf("%f\n", p->inf_fn_scale);
}

void set_seedRNG(int seed, params_t *p) {
    // random number seed
    time_t t;
    if (!seed) 
        seed = time(&t);
    p->stream = gsl_rng_alloc(gsl_rng_taus);    
    gsl_rng_set(p->stream, seed);
}

void output_result(FILE *fp, result_t *r, params_t *p, enum model_type type) {
    int y, m, month;
    char model_name[][20] = {"matrix", "deterministic", "stochastic", "parameters"};

    fprintf(fp, "%s max_age %d\n", model_name[type], p->max_age+1);
    
    if (type == Parameters) {
        fprintf(fp, "sim_years %d\n", p->sim_years);
        fprintf(fp, "flock_size %f\n", p->flock_size);
        fprintf(fp, "all_ewes_flock_size %f\n", p->all_ewes_flock_size);
        fprintf(fp, "max_age %d\n", p->max_age);
        fprintf(fp, "age_first_mating %d\n", p->age_first_mating);
        fprintf(fp, "pens %f\n", p->pens);
        fprintf(fp, "housing_period %d\n", p->housing_period);
        fprintf(fp, "beta_housed %f\n", p->beta_housed);
        fprintf(fp, "beta_field %f\n", p->beta_field);
        fprintf(fp, "L_mean %f\n", p->L_mean);
        fprintf(fp, "L_sd %f\n", p->L_sd);
        fprintf(fp, "S_mean %f\n", p->S_mean);
        fprintf(fp, "S_sd %f\n", p->S_sd);
        fprintf(fp, "inter_pen_contact_weight %f\n", p->inter_pen_contact_weight);
        fprintf(fp, "prob_dam_to_lamb %f\n", p->prob_dam_to_lamb);
        fprintf(fp, "%s\n", p->text);
    }
    else if (type == Next_Gen_Matrix) {
        fprintf(fp, "R0 %f\n", r->R0);
        fprintf(fp, "generation_time %f\n", r->generation_time);
        fprintf(fp, "doubling_time %f\n", r->doubling_time);
        fprintf(fp, "init_time %d\n", r->init_time);
    }
    else {
        for (y = 0; y < r->sim_years; y++)
            for (m = 0; m < 12; m++) {
                month = 12*y+m;
                fprintf(fp, "%d ", month);
                // for (i = 0; i <= p->max_age; i++) {
                //     fprintf(fp, "%.2f ", r->infected_by_age[month][i]/r->total_by_age[month][i]);
                //     fprintf(fp, "%.2f ", r->positive_by_age[month][i]/r->total_by_age[month][i]);
                //     fprintf(fp, "%.2f ", r->susceptible_by_age[month][i]/r->total_by_age[month][i]);
                //     // fprintf(fp, "%.2f ", r->infectiousness_by_age[month][i]);
                //     fprintf(fp, "%.2f ", r->foi_by_age[month][i]);
                //     fprintf(fp, "%.2f ", r->total_by_age[month][i]);
                //     fprintf(fp, "%.2f ", r->positive_by_age[month][i]);
                //     fprintf(fp, "%.2f ", r->infected_by_age[month][i]);
                //     fprintf(fp, "%.2f ", r->tsi[month][i]);
                // }
                // fprintf(fp, "%.2f ", r->positive[month]);
                // fprintf(fp, "%.2f ", r->infected[month]);
                // fprintf(fp, "%.2f ", r->resistant[month]);
                // fprintf(fp, "%.2f ", r->positive_flock[month]/r->flock_size[month]);
                // fprintf(fp, "%.2f ", r->positive_flock[month]);
                fprintf(fp, "%.2f ", r->infected_flock[month]/r->flock_size[month]);
                // fprintf(fp, "%.2f ", r->infected_flock[month]);
                // fprintf(fp, "%.2f ", r->seroconversions[month]);
                // fprintf(fp, "%.2f ", r->prop_pos[month]);
                // fprintf(fp, "%.2f ", r->mean_infectious_age[month]);
                fprintf(fp, "\n");
            }
        // for (y = 0; y < r->sim_years + p->max_age+1; y++)
        //     for (m = 0; m < p->max_age*12; m++)
        //         fprintf(fp, "%.2f ", r->infected_by_cohort[y][m]);
    }

}

double*** init_contact_weights(int max_age) {
    // int m, c;
    double ***contact_weights = doublematrix3(12, max_age+1, max_age+1);

    // initial weighting is set to zero. weights are set 
    // within age class (or cohort) weighting is 1 unless the cohort
    // mixes within pens with other cohorts 
    // for (m = April; m <= March; m++)
    //     for (c = 0; c <= max_age; c++)
    //         contact_weights[m][c][c] = 1;
    return contact_weights;
}

void print_contact_weights(double ***contact_weights, enum month m, params_t *p) {
    // called from ibm->Stochastic_Sim
    int i, j;
    int max_age;

    if (m == March) {
        printf("March housing\n");
        max_age = p->max_age-1;
    }
    else {
        printf("April pasture\n");
        max_age = p->max_age;
    }

    printf("  ");
    for (i = 0; i <= max_age; i++)
        printf("%d ", i);
    printf("\n");
        
    for (i = 0; i <= max_age; i++) {
        printf("%d ", i);
        for (j = 0; j <= max_age; j++)
            if (contact_weights[m][i][j] == 1)
                printf("* ");
            else if (m == March && contact_weights[m][i][j] > p->inter_pen_contact_weight)
                printf("o ");
            else if (contact_weights[m][i][j] > 0)
                printf(". ");
            else
                printf("  ");
        printf("\n");
    }
}

void run(params_t *p, params_t *init_p) {
    int i;
    double x;
    result_t *r;
    FILE *fp;
    char result_file[200] = "Results/";
    void *saved_cohorts;

    printf("\n######### %s\n", p->text);

    update_parameters(p);
    update_parameters(init_p);

    strcat(result_file, p->file);
    strcat(result_file, ".txt");
    fp = fopen(result_file, "w");

    output_result(fp, NULL, p, Parameters);

    // Stochastic simulation
#   pragma omp parallel for firstprivate(fp) private(x, saved_cohorts, r)
    for (i = 0; i < p->num_ibm_sims; i++) {
        do {
            saved_cohorts = NULL;
            r = Stochastic_Sim(init_p, &saved_cohorts);
            x = r->infected_flock[init_p->sim_years*12-1];
            result_free(r, Stochastic);
        } while(x < 150);
        r = Stochastic_Sim(p, &saved_cohorts);
#       pragma omp critical (output_result)
        output_result(fp, r, p, Stochastic);
        result_free(r, Stochastic);
    }    

    // Deterministic simulation
    r = Deterministic_Sim(p);
    output_result(fp, r, p, Deterministic);
    result_free(r, Deterministic);

    // Next generation matrix
    r = Next_Gen_Matrix_R0(p);
    printf("R0 = %.6f\n", r->R0);
    // printf("mean infected life expectancy is %.1f months\n", r->generation_time);
    output_result(fp, r, p, Next_Gen_Matrix);
    result_free(r, Next_Gen_Matrix);
    
    fclose(fp);
}

void run_single(params_t *p) {
    result_t *r;
    FILE *fp;
    char result_file[200] = "Results/";
    void *saved_cohorts;
    saved_cohorts = NULL; 

    printf("\n######### %s\n", p->text);

    update_parameters(p);

    strcat(result_file, p->file);
    strcat(result_file, ".txt");
    fp = fopen(result_file, "w");

    output_result(fp, NULL, p, Parameters);

    // Stochastic simulation
    r = Stochastic_Sim(p, &saved_cohorts);
    output_result(fp, r, p, Stochastic);
    result_free(r, Stochastic);

    // Deterministic simulation
    // r = Deterministic_Sim(p);
    // output_result(fp, r, p, Deterministic);
    // result_free(r, Deterministic);

    // Next generation matrix
    r = Next_Gen_Matrix_R0(p);
    printf("R0 = %.6f\n", r->R0);
    printf("mean infected life expectancy is %.1f months\n", r->generation_time);
    output_result(fp, r, p, Next_Gen_Matrix);
    result_free(r, Next_Gen_Matrix);
    
    fclose(fp);
}

void run_multiple(params_t *p) {
    result_t *r;
    FILE *fp;
    int i;
    char result_file[200] = "Results/";
    void *saved_cohorts;
    saved_cohorts = NULL; 

    printf("\n######### %s\n", p->text);

    update_parameters(p);

    strcat(result_file, p->file);
    strcat(result_file, ".txt");
    fp = fopen(result_file, "w");

    output_result(fp, NULL, p, Parameters);

    // Stochastic simulation
#   pragma omp parallel for firstprivate(fp) private(saved_cohorts, r)
    for (i = 0; i < p->num_ibm_sims; i++) {
        saved_cohorts = NULL; 
        r = Stochastic_Sim(p, &saved_cohorts);
#       pragma omp critical (output_result)
        output_result(fp, r, p, Stochastic);
        result_free(r, Stochastic);
    }    

    // r = Stochastic_Sim(p, &saved_cohorts);
    // output_result(fp, r, p, Stochastic);
    // result_free(r, Stochastic);

    // Deterministic simulation
    // r = Deterministic_Sim(p);
    // output_result(fp, r, p, Deterministic);
    // result_free(r, Deterministic);

    // Next generation matrix
    r = Next_Gen_Matrix_R0(p);
    printf("R0 = %.6f\n", r->R0);
    printf("mean infected life expectancy is %.1f months\n", r->generation_time);
    output_result(fp, r, p, Next_Gen_Matrix);
    result_free(r, Next_Gen_Matrix);
    
    fclose(fp);
}

void get_retain_probabilities(double *alpha, double *uninfected_cohort_size, double *lambs_per_dam, params_t *p) {
    int max_age = p->max_age;
    int age_first_mating = p->age_first_mating;
    int sale_month = p->sale_month;
    // double lambs_per_dam = p->lambs_per_dam;

    int i;
    double nlambs;
    double remove = uninfected_cohort_size[sale_month-1] - uninfected_cohort_size[sale_month];

    if (p->sale_strategy == Oldest) {
        // remove nlambs of oldest dams
        // number of dams surviving to this age time lambs_per_dam
        nlambs = uninfected_cohort_size[12*max_age + April] * lambs_per_dam[max_age];
        if (remove <= nlambs) {
            alpha[max_age] = 1. - remove/nlambs;
            remove = 0;
        }
        else {
            alpha[max_age] = 0;
            remove -= nlambs;
        }

        // remove equal proportions of nlambs from all other dam ages
        nlambs = 0;
        for (i = age_first_mating+1; i < max_age; i++)
            nlambs += uninfected_cohort_size[12*i + April] * lambs_per_dam[i];
        for (i = age_first_mating+1; i < max_age; i++)
            alpha[i] = 1. - remove / nlambs;
    }
    else
        for (i = age_first_mating+1; i <= max_age; i++)
            alpha[i] = uninfected_cohort_size[sale_month]/uninfected_cohort_size[sale_month-1];
}

double** litter_dist_by_age_data(int age_first_mating, int max_age) {
    // litter size (columns from 0 lambs to 4 lambs) by age (rows)
    // only twins born
    // double data[5][5] = {
    //     {0, 0, 1, 0, 0},
    //     {0, 0, 1, 0, 0},
    //     {0, 0, 1, 0, 0},
    //     {0, 0, 1, 0, 0},
    //     {0, 0, 1, 0, 0}
    // };
    // data from Andrew
    double data[5][5] = {
        {0.033457, 0.446097, 0.468401, 0.048327, 0.003717},
        {0.034043, 0.344681, 0.570213, 0.051064, 0.000000},
        {0.018182, 0.303030, 0.618182, 0.060606, 0.000000},
        {0.040404, 0.373737, 0.565657, 0.020202, 0.000000},
        {0.000000, 0.285714, 0.619048, 0.095238, 0.000000}
    };
    double **x = doublematrix(7, 5);
    int i, j;

    for (i = age_first_mating+1; i <= max_age; i++)
        for (j = 0; j < 5; j++)
            x[i][j] = data[i-2][j];
    return x;
}
#include "R0.h"

typedef struct Ewe {
    enum state state;
    enum sero_state sero_status;
    int tsi;
    int latent_period;
    int seroconversion_period;
    int dam_age_at_birth;
    double growth_rate;
    struct Ewe *dam;  // pointer to dam of this ewe
} ewe_t;

typedef struct {
    bool exists;
    int age_in_months;
    int num_ewes;
    ewe_t *ewes;
} cohort_t;

struct str {
    int value;
    int index;
};

static cohort_t* init_cohorts(int n_cohorts) {
    size_t c;
    cohort_t* cohorts = (cohort_t*) malloc(n_cohorts * sizeof(cohort_t));
    for (c = 0; c < n_cohorts; c++) {
        cohorts[c].exists = false;
        cohorts[c].age_in_months = 0;
    }
    return cohorts;
}

static cohort_t* get_cohort_of_age(int age_in_years, cohort_t *cohorts, int n_cohorts) {
    // get the cohort of age age_in_years
    size_t c;
    for (c = 0; c < n_cohorts; c++)
        if (cohorts[c].exists && cohorts[c].age_in_months / 12 == age_in_years)
            return &cohorts[c];
    return NULL;
}

static void cohorts_free(cohort_t *cohorts, int n_cohorts) {
    size_t c;
    for (c = 0; c < n_cohorts; c++)
        free(cohorts[c].ewes);
    free(cohorts);
}

static void set_cohort(cohort_t *cohort, int age_in_months, int num_ewes, double prop_resistant, gsl_rng *stream) {
    int i;
    cohort->exists = true;
    cohort->num_ewes = num_ewes;
    cohort->age_in_months = age_in_months;
    cohort->ewes = (ewe_t*) malloc(num_ewes * sizeof(ewe_t));
    for (i = 0; i < num_ewes; i++) {
        if (gsl_ran_bernoulli(stream, prop_resistant))
            cohort->ewes[i].state = Resistant;
        else
            cohort->ewes[i].state = Susceptible;
        cohort->ewes[i].sero_status = Negative;
    }
}

static int get_number_in_cohort_in_sero_state(enum sero_state s, cohort_t *cohort) {
    // get number of ewes in cohort in sero_state s
    int i, count = 0;
    for (i = 0; i < cohort->num_ewes; i++)
        if (cohort->ewes[i].state != Removed && cohort->ewes[i].sero_status == s)
        count++;
    return count;
}

static int get_number_in_cohort_in_state(enum state s, cohort_t *cohort) {
    // get number of ewes in cohort in state s
    if (s == Any)
        return cohort->num_ewes - get_number_in_cohort_in_state(Removed, cohort);
    
    if (s == Infected)
        return get_number_in_cohort_in_state(Exposed, cohort) + get_number_in_cohort_in_state(Infectious, cohort);

    // otherwise get number in a single state, Susceptible, Exposed, Infectious, Resistant
    int i, count = 0;
    for (i = 0; i < cohort->num_ewes; i++)
        if (cohort->ewes[i].state == s)
            count++;
    return count;
}

static int get_number_infected_for(int y, cohort_t *cohorts, int n_cohorts, int max_age) {
    // get number of ewes infected for i years
    int i, n, count = 0;
    cohort_t *cohort;

    for (i = 0; i <= max_age; i++) {
        cohort = get_cohort_of_age(i, cohorts, n_cohorts);
        if (cohort)
            for (n = 0; n < cohort->num_ewes; n++)
                if (cohort->ewes[n].state == Exposed || cohort->ewes[n].state == Infectious)
                    if (cohort->ewes[n].tsi/12 == y)
                        count++;
    }
    return count;
}

static int get_number_in_sero_state(enum sero_state s, cohort_t *cohorts, int n_cohorts, int max_age) {
    // get number of ewes in sero_state s
    int i, count = 0;
    cohort_t *cohort;

    for (i = 0; i <= max_age; i++) {
        cohort = get_cohort_of_age(i, cohorts, n_cohorts);
        if (cohort)
            count += get_number_in_cohort_in_sero_state(s, cohort);
    }
    return count;
}

static int get_number_in_state(enum state s, cohort_t *cohorts, int n_cohorts, int max_age) {
    // get number of ewes in state s
    int i, count = 0;
    cohort_t *cohort;

    for (i = 0; i <= max_age; i++) {
        cohort = get_cohort_of_age(i, cohorts, n_cohorts);
        if (cohort)
            count += get_number_in_cohort_in_state(s, cohort);
    }
    return count;
}

static int get_number_mated_in_sero_state(enum sero_state s, cohort_t *cohorts, int n_cohorts, int age_first_mating, int max_age) {
    // count the number of ewes that will be, or have been, mated this year
    // ie do not count the oldest ewes which will be sold before being mated
    int i, count = 0;
    cohort_t *cohort;
    
    for (i = age_first_mating; i < max_age; i++) {
        cohort = get_cohort_of_age(i, cohorts, n_cohorts);
        if (cohort)
            count += get_number_in_cohort_in_sero_state(s, cohort);
    }
    return count;
}

static int get_number_mated_in_state(enum state s, cohort_t *cohorts, int n_cohorts, int age_first_mating, int max_age) {
    // count the number of ewes that will be, or have been, mated this year
    // ie do not count the oldest ewes which will be sold before being mated
    int i, count = 0;
    cohort_t *cohort;
    
    for (i = age_first_mating; i < max_age; i++) {
        cohort = get_cohort_of_age(i, cohorts, n_cohorts);
        if (cohort)
            count += get_number_in_cohort_in_state(s, cohort);
    }
    return count;
}

static double get_cohort_infectiousness(cohort_t *cohort, params_t *p) {
    // get infectiousness of all infectious ewes in a cohort
    if (cohort == NULL)
        return 0;
    
    int i;
    ewe_t *ewe;
    double infectiousness = 0;

    for (i = 0; i < cohort->num_ewes; i++) {
        ewe = &cohort->ewes[i]; 
        if (ewe->state == Infectious)
            infectiousness += p->inf_fn_scale * (*p->inf_fn)(ewe->tsi, p);
    }
    return infectiousness;
}

static void infect_ewes_in_cohort(int n, cohort_t *cohort, params_t *p, enum state s) {
    // randomly infect n Susceptible ewes in this cohort
    // ewes can be made immediately infectious (ie first ewe) by passing s=Infectious
    int i;
    ewe_t *ewe;

    while (n > 0) {
        // choose a random ewe in this cohort
        i = gsl_rng_uniform_int(p->stream, cohort->num_ewes);
        ewe = &cohort->ewes[i]; 
        if (ewe->state == Susceptible) {
            if (s == Exposed) {
                ewe->state = Exposed;
                ewe->tsi = 0;
                ewe->latent_period = lrint(gsl_ran_gamma(p->stream, p->L_a, p->L_b));
                ewe->seroconversion_period = lrint(gsl_ran_gamma(p->stream, p->S_a, p->S_b));
            }
            else {
                // ewe is infectious so set time since infection to mean latent period
                ewe->state = Infectious;
                ewe->tsi = lrint(p->L_mean);
                ewe->latent_period = lrint(p->L_mean);
                ewe->seroconversion_period = lrint(p->S_mean);
            }
            n--;
        }
    }
}

static void remove_ewes_from_cohort(int n, cohort_t *cohort, gsl_rng *stream) {
    // randomly remove existing ewes from this cohort
    int i;
    ewe_t *ewe;

    while (n > 0) {
        i = gsl_rng_uniform_int(stream, cohort->num_ewes);
        ewe = &cohort->ewes[i]; 
        if (ewe->state != Removed) {
            ewe->state = Removed;
            n--;
        }
    }
}

static int cmp(const void *a, const void *b)
{
    struct str *a1 = (struct str *)a;
    struct str *a2 = (struct str *)b;
    if ((*a1).value > (*a2).value)
        return 1;
    else if ((*a1).value < (*a2).value)
        return -1;
    else
        return 0;
}

static int sell_lambs(int n, cohort_t *lamb_cohort, params_t *p) {
    int i;
    int positive = 0; // the number of lambs sold with seropositive dams
    int nlambs = lamb_cohort->num_ewes;
    ewe_t *lamb, *dam;

    if (p->sale_strategy == Poor) {
        // preferentially remove lambs of dams with the smallest growth rate
        struct str *objects = (struct str*) malloc(nlambs * sizeof(struct str));
        for (i = 0; i < nlambs; i++) {
            objects[i].value = lamb_cohort->ewes[i].growth_rate;
            objects[i].index = i;
        }
        // sort lambs with ascending growth rate
        qsort(objects, nlambs, sizeof(objects[0]), cmp);

        // sell the n lambs with the smallest growth rate
        for (n--; n >= 0; n--) {
            lamb = &lamb_cohort->ewes[objects[n].index];
            lamb->state = Removed;
            dam = lamb->dam;
            if (dam != NULL && dam->sero_status == Positive)
                positive++; 
        }
        free(objects);
    }
    else if (p->sale_strategy == Oldest) {
        // first remove lambs of oldest dams
        for (i = 0; i < nlambs; i++) {
            lamb = &lamb_cohort->ewes[i];
            if (lamb->state != Removed && lamb->dam_age_at_birth == p->max_age) { 
                lamb->state = Removed;
                dam = lamb->dam;
                if (dam != NULL && dam->sero_status == Positive)
                    positive++; 
                n--;
            }
        }
    }

    // randomly remove the remaining lambs of dams of all other ages
    remove_ewes_from_cohort(n, lamb_cohort, p->stream);
    return positive;
}

static int age_infected_ewes(cohort_t *cohort) {
    // age all ewes in this cohort and set them Infectious if time since infection
    // equals latent period
    int i;
    int count_sero = 0;
    ewe_t *ewe;

    for (i = 0; i < cohort->num_ewes; i++) {
        ewe = &cohort->ewes[i]; 
        // if reached latent period ewe becomes infectious
        if (ewe->state == Exposed)
            if (ewe->tsi >= ewe->latent_period)
                ewe->state = Infectious;

        if (ewe->state == Exposed || ewe->state == Infectious) {
            // if infected ewe reaches seroconversion period ewe becomes sero positive
            if (ewe->sero_status == Negative && ewe->tsi >= ewe->seroconversion_period) {
                ewe->sero_status = Positive;
                count_sero++;
            }

            // increment infected ewe time since infection
            ewe->tsi++;
        }
    }
    // return number of seroconverted ewes this month
    return count_sero;
}

static void add_lamb(ewe_t *dam, int dam_age, cohort_t *lamb_cohort, enum state s, params_t *p) {
    lamb_cohort->num_ewes += 1;
    size_t n = lamb_cohort->num_ewes;
    lamb_cohort->ewes = (ewe_t*) realloc(lamb_cohort->ewes, n * sizeof(ewe_t));
    ewe_t *lamb = &lamb_cohort->ewes[n-1];
    lamb->state = s;
    lamb->sero_status = Negative;
    lamb->dam_age_at_birth = dam_age;
    lamb->dam = dam;
    lamb->tsi = 0;
    lamb->growth_rate = p->gr_tsi_grad * dam->tsi + gsl_ran_gaussian(p->stream, 1.0);
    if (s == Exposed) {
        lamb->latent_period = gsl_ran_gamma(p->stream, p->L_a, p->L_b);
        lamb->seroconversion_period = gsl_ran_gamma(p->stream, p->S_a, p->S_b);
    }
}

static void load_cohorts(cohort_t *saved_cohorts, cohort_t *cohorts, int n_cohorts, int max_age) {
    int i, j, c;
    cohort_t *cohort;

    for (i = 0; i <= max_age; i++) {
        c = max_age - i;
        cohort = &cohorts[c];
        cohort->exists = saved_cohorts[i].exists;
        cohort->age_in_months = saved_cohorts[i].age_in_months;
        cohort->num_ewes = saved_cohorts[i].num_ewes;
        cohort->ewes = (ewe_t*) malloc(saved_cohorts[i].num_ewes * sizeof(ewe_t));
        for (j = 0; j < saved_cohorts[i].num_ewes; j++) {
            cohort->ewes[j].state = saved_cohorts[i].ewes[j].state;
            cohort->ewes[j].sero_status = saved_cohorts[i].ewes[j].sero_status;
            cohort->ewes[j].tsi = saved_cohorts[i].ewes[j].tsi;
            cohort->ewes[j].latent_period = saved_cohorts[i].ewes[j].latent_period;
            cohort->ewes[j].seroconversion_period = saved_cohorts[i].ewes[j].seroconversion_period;
            cohort->ewes[j].dam_age_at_birth = saved_cohorts[i].ewes[j].dam_age_at_birth;
            cohort->ewes[j].growth_rate = saved_cohorts[i].ewes[j].growth_rate;
            cohort->ewes[j].dam = NULL;
        }
    }
}

static cohort_t *save_cohorts(cohort_t *cohorts, int n_cohorts, int max_age) {
    int i, j;
    cohort_t *cohort;
    cohort_t *saved_cohorts = init_cohorts(max_age+1);
    
    for (i = 0; i <= max_age; i++) {
        cohort = get_cohort_of_age(i, cohorts, n_cohorts);
        if (cohort) {
            saved_cohorts[i].exists = cohort->exists;
            saved_cohorts[i].age_in_months = cohort->age_in_months;
            saved_cohorts[i].num_ewes = cohort->num_ewes;
            saved_cohorts[i].ewes = (ewe_t*) malloc(cohort->num_ewes * sizeof(ewe_t));
            for (j = 0; j < cohort->num_ewes; j++) {
                saved_cohorts[i].ewes[j].state = cohort->ewes[j].state;
                saved_cohorts[i].ewes[j].sero_status = cohort->ewes[j].sero_status;
                saved_cohorts[i].ewes[j].tsi = cohort->ewes[j].tsi;
                saved_cohorts[i].ewes[j].latent_period = cohort->ewes[j].latent_period;
                saved_cohorts[i].ewes[j].seroconversion_period = cohort->ewes[j].seroconversion_period;
                saved_cohorts[i].ewes[j].dam_age_at_birth = cohort->ewes[j].dam_age_at_birth;
                saved_cohorts[i].ewes[j].growth_rate = cohort->ewes[j].growth_rate;
                saved_cohorts[i].ewes[j].dam = NULL;
            }
        }
    }
    return saved_cohorts;
}

result_t* Stochastic_Sim(params_t *p, void **saved_cohorts) {
    // Stochastic simulation of an outbreak starting with a single infected lamb 
    int max_age = p->max_age;
    int sim_years = p->sim_years;
    int age_first_mating = p->age_first_mating;
    enum month sale_month = p->sale_month;
    enum month mating_month = p->mating_month;
    // enum month weaning_month = p->weaning_month;
    double prob_dam_to_lamb = p->prob_dam_to_lamb;
    // double beta_housed = p->beta_housed;
    // double L_a = p->L_a;
    // double L_b = p->L_b;
    gsl_rng *stream = p->stream;

    double *p_infectious = prob_infectious(p);
    double *removal_rate = natural_deaths(p->n_age_classes);
    double *infection_rate = monthly_infection_rate(p->housing_period, p->beta_housed, p->beta_field);
    double *uninfected_cohort_size = cohort_size(p, removal_rate);
    double ***contact_weights = (*p->contact_weight_func)(p);

    // ewe_t *lamb;
    ewe_t *dam;
    enum month m;
    int y, c, i, a, n, seroconversions;
    int age_in_months, month;
    int oldest_cohort, youngest_cohort;
    double prop_pos = 0;
    double p_alive, p_infection, foi;
    double cohort_infectiousness[max_age+1];
    double retained = uninfected_cohort_size[sale_month];

    int n_cohorts = sim_years+max_age+1;
    cohort_t *cohorts = init_cohorts(n_cohorts);
    cohort_t *cohort, *dam_cohort, *lamb_cohort;
    result_t *r = init_results(sim_years, max_age, 0);

    double **litter_dist_by_age = litter_dist_by_age_data(age_first_mating, max_age);
    gsl_ran_discrete_t *g[max_age+1];
    for (a = age_first_mating+1; a <= max_age; a++)
        g[a] = gsl_ran_discrete_preproc(5, litter_dist_by_age[a]);

    // print_contact_weights(contact_weights, March, p);

    // set cohort ages, cohort 0 is the oldest, 
    // cohort max_age is the youngest at age 0
    oldest_cohort = 0;
    youngest_cohort = max_age;

    // set up all cohorts at start of new year (start of April) with none infected
    if (*saved_cohorts == NULL)
        for (c = oldest_cohort; c <= youngest_cohort; c++) {
            a = 12*(max_age-c);
            set_cohort(&cohorts[c], a, uninfected_cohort_size[a], p->prop_resistant, stream);
        }
    else
        load_cohorts(*saved_cohorts, cohorts, n_cohorts, max_age);

    // simulate over each year starting in April
    for (y = 0; y < sim_years; y++) {
        // set number of seroconversions each year at start of year to 0
        seroconversions = 0;

        for (m = April; m < March; m++) {
            month = 12*y+m;

            // in the first mating month make a 1 year old in mated flock infectious
            if (*saved_cohorts == NULL && month == mating_month)
                infect_ewes_in_cohort(1, &cohorts[youngest_cohort-1], p, Infectious);

            // calculate number of infectious ewes and other stats at the start of this month
            r->positive[month] = get_number_in_sero_state(Positive, cohorts, n_cohorts, max_age);
            r->susceptible[month] = get_number_in_state(Susceptible, cohorts, n_cohorts, max_age);
            r->infected[month] = get_number_in_state(Infected, cohorts, n_cohorts, max_age);
            r->infectious[month] = get_number_in_state(Infectious, cohorts, n_cohorts, max_age);
            r->resistant[month] = get_number_in_state(Resistant, cohorts, n_cohorts, max_age);
            r->total[month] = get_number_in_state(Any, cohorts, n_cohorts, max_age);
            r->positive_flock[month] = get_number_mated_in_sero_state(Positive, cohorts, n_cohorts, age_first_mating, max_age);
            r->infected_flock[month] = get_number_mated_in_state(Infected, cohorts, n_cohorts, age_first_mating, max_age);
            r->infectious_flock[month] = get_number_mated_in_state(Infectious, cohorts, n_cohorts, age_first_mating, max_age);
            r->flock_size[month] = get_number_mated_in_state(Any, cohorts, n_cohorts, age_first_mating, max_age);
            r->seroconversions[month] = seroconversions;
            r->prop_pos[month] = prop_pos;

            for (i = 0; i <= max_age; i++) {
                cohort = get_cohort_of_age(i, cohorts, n_cohorts);
                if (cohort) {
                    r->positive_by_age[month][i] = get_number_in_cohort_in_sero_state(Positive, cohort);
                    r->susceptible_by_age[month][i] = get_number_in_cohort_in_state(Susceptible, cohort);
                    r->infected_by_age[month][i] = get_number_in_cohort_in_state(Infected, cohort);
                    r->infectious_by_age[month][i] = get_number_in_cohort_in_state(Infectious, cohort);
                    r->resistant_by_age[month][i] = get_number_in_cohort_in_state(Resistant, cohort);
                    r->total_by_age[month][i] = get_number_in_cohort_in_state(Any, cohort);
                    r->tsi[month][i] = get_number_infected_for(i, cohorts, n_cohorts, max_age);
                } else {
                    r->positive_by_age[month][i] = 0;
                    r->susceptible_by_age[month][i] = 0;
                    r->infected_by_age[month][i] = 0;
                    r->infectious_by_age[month][i] = 0;
                    r->resistant_by_age[month][i] = 0;
                    r->total_by_age[month][i] = 0;
                    r->tsi[month][i] = 0;
                }
            }
            
            for (c = oldest_cohort; c <= youngest_cohort; c++) {
                cohort = &cohorts[c];
                age_in_months = cohort->age_in_months;
                r->infected_by_cohort[c][age_in_months] = get_number_in_cohort_in_state(Infected, cohort);
            }

            // get infectiousness of each cohort
            for (i = 0; i <= max_age; i++) {
                cohort = get_cohort_of_age(i, cohorts, n_cohorts);
                cohort_infectiousness[i] = get_cohort_infectiousness(cohort, p);
                r->infectiousness_by_age[month][i] = cohort_infectiousness[i];
            }

            // infections during this month for each cohort
            for (c = oldest_cohort; c <= youngest_cohort; c++) {
                cohort = &cohorts[c];
                age_in_months = cohort->age_in_months;
                // prob a ewe is still alive at end of month
                p_alive = exp(-removal_rate[age_in_months]);

                // force of infection experienced by this cohort  
                foi = 0;
                for (i = 0; i <= max_age; i++)
                    foi += contact_weights[m][i][age_in_months/12] * cohort_infectiousness[i];
                foi *= infection_rate[m];
                r->foi_by_age[month][age_in_months/12] = foi;
                p_infection = 1-exp(-foi);

                // natural death of random ewes
                n = get_number_in_cohort_in_state(Any, cohort);
                n = gsl_ran_binomial(stream, 1-p_alive, n);
                remove_ewes_from_cohort(n, cohort, stream);

                // number of newly infected
                n = get_number_in_cohort_in_state(Susceptible, cohort);
                // printf("%d %d %d %d ", y, m, c, n);
                n = gsl_ran_binomial(stream, p_infection, n);
                // printf("%d ", n);
                infect_ewes_in_cohort(n, cohort, p, Exposed);
                n = get_number_in_cohort_in_state(Susceptible, cohort);
                // printf("%d %f\n", n, p_infection);

                // ignore this because we assume that all modes of virus transmission from dam 
                // to lamb is modelled via prob_dam_to_lamb
                // lamb infection from dam before weaning
                // maternal transmission via exhalation of virus during weaning
                // if (y && c == youngest_cohort && m < weaning_month) {
                //     for (i = 0; i < cohort->num_ewes; i++) {
                //         lamb = &cohort->ewes[i];
                //         dam = lamb->dam;
                //         if (lamb->state == Susceptible && dam->state == Infectious) {
                //             foi = p->inf_fn_scale * (*p->inf_fn)(dam->tsi, p) * p->beta_housed;
                //             if (gsl_ran_bernoulli(stream, 1-exp(-foi))) {
                //                 lamb->state = Exposed;
                //                 lamb->tsi = 0;
                //                 lamb->latent_period = lrint(gsl_ran_gamma(stream, L_a, L_b));
                //             }
                //         }
                //     }
                // }
                // increment age for start of next month
                seroconversions += age_infected_ewes(cohort);
                cohort->age_in_months++;
            }

            if (m == sale_month-1) {
                // sell enough lambs to maintain a stable flock size
                lamb_cohort = &cohorts[youngest_cohort];
                n = get_number_in_cohort_in_state(Any, lamb_cohort) - lrint(retained);
                prop_pos = (float) sell_lambs(n, lamb_cohort, p) / (float) n;
                // print()

                // sell ewes with high tsi if sale strategy is Poor
                // if (p->sale_strategy == Poor)
                //     for (c = oldest_cohort; c <= youngest_cohort; c++) {
                //         cohort = &cohorts[c];
                //         for (i = 0; i < cohort->num_ewes; i++) {
                //             dam = &cohort->ewes[i];
                //             if (dam->state == Exposed && dam->tsi > p->tsi_removal) {
                //                 dam->state = Removed;
                //             }
                //         }
                //     }

                // sell oldest cohort
                cohorts[oldest_cohort].exists = false;
                oldest_cohort++;
            }
        }

        // Now at beginning of April, add a new cohort of lambs
        youngest_cohort++;
        lamb_cohort = &cohorts[youngest_cohort];
        lamb_cohort->exists = true;
        lamb_cohort->num_ewes = 0;
        lamb_cohort->age_in_months = 0;
        lamb_cohort->ewes = NULL;

        // loop through all dams and create new lambs
        for (a = age_first_mating+1; a <= max_age; a++) {
            // get cohort of dams
            dam_cohort = get_cohort_of_age(a, cohorts, n_cohorts);
            // loop through all existing dams
            for (n = 0; n < dam_cohort->num_ewes; n++) {
                dam = &dam_cohort->ewes[n];
                if (dam->state != Removed) {
                    for (i = 0; i < gsl_ran_discrete(stream, g[a]); i++)
                        // if dam is resitant so is lamb
                        if (dam->state == Resistant)
                            add_lamb(dam, a, lamb_cohort, Resistant, p);
                        else if (dam->state == Infectious && gsl_ran_bernoulli(stream, prob_dam_to_lamb))
                            // maternal infection via infected colostrum
                            add_lamb(dam, a, lamb_cohort, Exposed, p);
                        else
                            add_lamb(dam, a, lamb_cohort, Susceptible, p);
                }
            }
        }
    }

    // printf("\n ##### Stochastic simulation #####\n");
    // for (y = 0; y < sim_years; y++) {
    //     for (m = 0; m < 12; m++) {
    //         month = 12*y+m;
    //         if (1||m == mating_month) {
    //             printf("% 3d % 3d | ", y, m);
    //             if (m == sale_month)
    //                 printf("S | ");
    //             else if (m == mating_month)
    //                 printf("M | ");
    //             else if (m == 0)
    //                 printf("B | ");
    //             else
    //                 printf("  | ");
    //             for (i = 0; i <= max_age; i++)
    //                 printf("% 4.0f ", r->susceptible_by_age[month][i]);
    //             printf("= ");
    //             printf("% 4.0f | ", r->susceptible[month]);
    //             for (i = 0; i <= max_age; i++)
    //                 printf("% 3.0f ", r->infected_by_age[month][i]);
    //             printf("= ");
    //             printf("% 4.0f | ", r->infected[month]);
    //             for (i = 0; i <= max_age; i++)
    //                 printf("% 3.0f ", r->infectious_by_age[month][i]);
    //             printf("= ");
    //             printf("% 4.0f | ", r->infectious[month]);
    //             for (i = 0; i <= max_age; i++)
    //                 printf("% 3.0f ", r->total_by_age[month][i]);
    //             printf("= ");
    //             printf("% 4.0f | ", r->total[month]);
    //             printf("% 3.0f\n", r->flock_size[month]);
    //         }
    //     }
    // }
    // exit(0);
    
    if (*saved_cohorts == NULL)
        *saved_cohorts = save_cohorts(cohorts, n_cohorts, max_age);
    else
        cohorts_free(*saved_cohorts, max_age+1);

    // Tidy up
    free(litter_dist_by_age[0]);
    free(litter_dist_by_age);
    free(p_infectious);
    free(removal_rate);
    free(infection_rate);
    free(uninfected_cohort_size);
    free(contact_weights[0][0]);
    free(contact_weights[0]);
    free(contact_weights);
    cohorts_free(cohorts, n_cohorts);
    for (a = age_first_mating+1; a <= max_age; a++)
        gsl_ran_discrete_free(g[a]);

    return r;
}

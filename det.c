#include "R0.h"

typedef struct {
    bool exists;
    int age_in_months;
    double S;
    double *I;
    double *S_lamb;
    double **I_lamb;
} cohort_t;

static cohort_t* init_cohorts(int n_cohorts, int n_age_classes, int max_age) {
    size_t c;
    cohort_t *cohorts = (cohort_t*) malloc(n_cohorts * sizeof(cohort_t));
    for (c = 0; c < n_cohorts; c++) {
        cohorts[c].exists = false;
        cohorts[c].age_in_months = 0;
        cohorts[c].S = 0;
        cohorts[c].I = doublevector(n_age_classes);
        cohorts[c].S_lamb = doublevector(max_age+1);
        cohorts[c].I_lamb = doublematrix(max_age+1, n_age_classes);
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
    for (c = 0; c < n_cohorts; c++) {
        free(cohorts[c].I);
        free(cohorts[c].S_lamb);
        free(cohorts[c].I_lamb[0]);
        free(cohorts[c].I_lamb);
    }
    free(cohorts);
}

static void set_cohort(cohort_t *cohort, int age_in_months, double num_ewes) {
    int i;
    cohort->exists = true;
    cohort->age_in_months = age_in_months;
    cohort->S = num_ewes;
    for (i = 0; i <= age_in_months; i++)
        cohort->I[i] = 0;
}

static double get_number_in_cohort_in_state(enum state s, cohort_t *cohort, double *p_infectious) {
    // get number of ewes in cohort in state s
    int a;
    double count = 0;

    if (s == Susceptible || s == Any)
        count = cohort->S;
        
    for (a = 0; a <= cohort->age_in_months; a++) {
        if (s == Exposed)
            count += cohort->I[a] * (1-p_infectious[a]);
            
        else if (s == Infectious)
            count += cohort->I[a] * p_infectious[a];

        else if (s == Infected || s == Any)
            count += cohort->I[a];
    }
    return count;
}

static double get_number_in_state(enum state s, cohort_t *cohorts, int n_cohorts, int max_age, double *p_infectious) {
    int i;
    double count = 0;
    cohort_t *cohort;

    for (i = 0; i <= max_age; i++) {
        cohort = get_cohort_of_age(i, cohorts, n_cohorts);
        if (cohort)
            count += get_number_in_cohort_in_state(s, cohort, p_infectious);
    }
    return count;
}

static double get_number_mated_in_state(enum state s, cohort_t *cohorts, int n_cohorts, int age_first_mating, int max_age, double *p_infectious) {
    // count the number of ewes that will be, or have been, mated this year
    // ie do not count the oldest ewes which will be sold before being mated
    int i;
    double count = 0;
    cohort_t *cohort;

    for (i = age_first_mating; i < max_age; i++) {
        cohort = get_cohort_of_age(i, cohorts, n_cohorts);
        if (cohort)
            count += get_number_in_cohort_in_state(s, cohort, p_infectious);
    }
    return count;
}

static double get_cohort_infectiousness(cohort_t *cohort, double *p_infectious, params_t *p) {
    // get infectiousness of cohort: number of infectious ewes * infectiousness since infection
    if (cohort == NULL)
        return 0;
        
    size_t tsi;
    double infectiousness = 0;
    
    for (tsi = 0; tsi <= cohort->age_in_months; tsi++)
        infectiousness += cohort->I[tsi] * p_infectious[tsi] * p->inf_fn_scale * (*p->inf_fn)(tsi, p);
    return infectiousness;
}

static double mean_infectious_age(cohort_t *cohorts, int n_cohorts, int max_age, double *p_infectious) {
    int i;
    double n, count = 0, mean_age = 0;
    cohort_t *cohort;

    for (i = 0; i <= max_age; i++) {
        cohort = get_cohort_of_age(i, cohorts, n_cohorts);
        if (cohort) {
            n = get_number_in_cohort_in_state(Infectious, cohort, p_infectious);
            count += n;
            mean_age += n * cohort->age_in_months;
        }
    }
    return mean_age/count;
}


result_t* Deterministic_Sim(params_t *p) {
    // Deterministic simulation of an outbreak starting with a single infectious gimmer
    int max_age = p->max_age;
    int sim_years = p->sim_years;
    int age_first_mating = p->age_first_mating;
    int n_age_classes = p->n_age_classes;
    enum month sale_month = p->sale_month;
    enum month mating_month = p->mating_month;
    enum month weaning_month = p->weaning_month;
    double prob_dam_to_lamb = p->prob_dam_to_lamb;

    double *p_infectious = prob_infectious(p);
    double *removal_rate = natural_deaths(p->n_age_classes);
    double *infection_rate = monthly_infection_rate(p->housing_period, p->beta_housed, p->beta_field);
    double *uninfected_cohort_size = cohort_size(p, removal_rate);
    double ***contact_weights = (*p->contact_weight_func)(p);

    enum month m;
    int y, c, i, a, da;
    int age_in_months, month, oldest_cohort, youngest_cohort;
    double alpha[max_age+1];
    double lambs_per_dam[max_age+1];
    double foi, p_alive, p_infected, S0, lamb_foi;
    double dams, infectious_dams;
    double cohort_infectiousness[max_age+1];
    double *Lambda = cumulative_removal_rate(n_age_classes, removal_rate);
    double **litter_dist_by_age = litter_dist_by_age_data(age_first_mating, max_age);

    int n_cohorts = sim_years+max_age+1;
    cohort_t *cohorts = init_cohorts(n_cohorts, n_age_classes, max_age);
    cohort_t *cohort, *lamb_cohort, *dam_cohort;
    result_t *r = init_results(sim_years, max_age, 0);

    // setup up age-specific lambs born per dam
    for (a = age_first_mating+1; a <= max_age; a++) {
        lambs_per_dam[a] = 0;
        for (i = 0; i < 5; i++)
            lambs_per_dam[a] += i * litter_dist_by_age[a][i];
    };

    get_retain_probabilities(alpha, uninfected_cohort_size, lambs_per_dam, p);

    // set cohort ages, cohort 0 is the oldest, 
    // cohort max_age is the youngest at age 0
    oldest_cohort = 0;
    youngest_cohort = max_age;

    // set up all cohorts at start of new year (start of April) with none infected
    for (c = oldest_cohort; c <= youngest_cohort; c++) {
        a = 12*(max_age-c);
        set_cohort(&cohorts[c], a, uninfected_cohort_size[a]);
    }

    // setup initial lambs by dam age
    lamb_cohort = &cohorts[youngest_cohort];
    for (i = age_first_mating+1; i <= max_age; i++) {
        dam_cohort = get_cohort_of_age(i, cohorts, n_cohorts);
        lamb_cohort->S_lamb[i] = get_number_in_cohort_in_state(Any, dam_cohort, NULL) * lambs_per_dam[i];
        lamb_cohort->I_lamb[i][0] = 0;
    }


    // simulate over each year starting in April
    for (y = 0; y < sim_years; y++) {
        // initialise the cumulative force of infection for this year
        for (m = April; m <= March; m++) {
            month = 12*y+m;

            // in the first mating month introduce a just-infectious gimmer
            if (month == mating_month) {
                cohorts[youngest_cohort-1].S -= 1;
                cohorts[youngest_cohort-1].I[12+month] = 1;
            }

            // calculate number of infectious ewes and other stats at the start of this month
            r->mean_infectious_age[month] = mean_infectious_age(cohorts, n_cohorts, max_age, p_infectious);
            r->susceptible[month] = get_number_in_state(Susceptible, cohorts, n_cohorts, max_age, NULL);
            r->infected[month] = get_number_in_state(Infected, cohorts, n_cohorts, max_age, NULL);
            r->infectious[month] = get_number_in_state(Infectious, cohorts, n_cohorts, max_age, p_infectious);
            r->total[month] = get_number_in_state(Any, cohorts, n_cohorts, max_age, NULL);
            r->infected_flock[month] = get_number_mated_in_state(Infected, cohorts, n_cohorts, age_first_mating, max_age, NULL);
            r->infectious_flock[month] = get_number_mated_in_state(Infectious, cohorts, n_cohorts, age_first_mating, max_age, p_infectious);
            r->flock_size[month] = get_number_mated_in_state(Any, cohorts, n_cohorts, age_first_mating, max_age, NULL);
            for (i = 0; i <= max_age; i++) {
                cohort = get_cohort_of_age(i, cohorts, n_cohorts);
                if (cohort) {
                    r->susceptible_by_age[month][i] = get_number_in_cohort_in_state(Susceptible, cohort, NULL);
                    r->infected_by_age[month][i] = get_number_in_cohort_in_state(Infected, cohort, NULL);
                    r->infectious_by_age[month][i] = get_number_in_cohort_in_state(Infectious, cohort, p_infectious);
                    r->total_by_age[month][i] = get_number_in_cohort_in_state(Any, cohort, NULL);
                } else {
                    r->susceptible_by_age[month][i] = 0;
                    r->infected_by_age[month][i] = 0;
                    r->infectious_by_age[month][i] = 0;
                    r->total_by_age[month][i] = 0;
                }
            }

            // get cumulative infectiousness of each cohort
            for (i = 0; i <= max_age; i++) {
                cohort = get_cohort_of_age(i, cohorts, n_cohorts);
                cohort_infectiousness[i] = get_cohort_infectiousness(cohort, p_infectious, p);
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

                // lamb cohort treated separately due to age of dam
                if (c == youngest_cohort) {
                    // set the totals to zero
                    cohort->S = 0;
                    for (a = 0; a <= age_in_months; a++)
                        cohort->I[a] = 0;
                        
                    for (da = age_first_mating+1; da <= max_age; da++) {
                        lamb_foi = foi;
                        // additional force of infection from infectious dams before weaning
                        if (y && m < weaning_month)
                            lamb_foi += p->beta_housed * cohort_infectiousness[da] / r->total_by_age[month][da];
                        r->foi_by_age[month][age_in_months/12] = lamb_foi;

                        p_infected = 1-exp(-lamb_foi);
                        S0 = cohort->S_lamb[da];
                        // natural death and infection of susceptibles
                        cohort->S_lamb[da] = S0 * p_alive * (1 - p_infected);
                        // natural death and ageing of infecteds

                        for (a = min(n_age_classes-1, age_in_months+1); a > 0; a--)
                            cohort->I_lamb[da][a] = cohort->I_lamb[da][a-1] * p_alive;
                        // new infecteds including natural death
                        cohort->I_lamb[da][0] = S0 * p_alive * p_infected;

                        // add to the totals
                        cohort->S += cohort->S_lamb[da];
                        for (a = 0; a <= min(n_age_classes-1, age_in_months+1); a++)
                            cohort->I[a] += cohort->I_lamb[da][a];
                    }
                }
                else {
                    // otherwise this cohort comprises ewes
                    p_infected = 1-exp(-foi);
                    S0 = cohort->S;
                    // natural death and infection of susceptibles
                    cohort->S = S0 * p_alive * (1 - p_infected);
                    // natural death and ageing of infecteds
                    for (a = min(n_age_classes-1, age_in_months+1); a > 0; a--)
                        cohort->I[a] = cohort->I[a-1] * p_alive;
                    // new infecteds including natural death
                    cohort->I[0] = S0 * p_alive * p_infected;
                }

                // increment age for start of next month
                cohort->age_in_months++;
            }

            // sale of lambs and oldest ewes for the start of the sale month
            if (m == sale_month-1) {
                cohort = &cohorts[youngest_cohort];
                // set the totals to zero
                cohort->S = 0;
                for (a = 0; a <= sale_month; a++)
                    cohort->I[a] = 0;

                for (da = age_first_mating+1; da <= max_age; da++) {
                    cohort->S_lamb[da] *= alpha[da];
                    for (a = 0; a <= sale_month; a++)
                        cohort->I_lamb[da][a] *= alpha[da];

                    // add to the totals
                    cohort->S += cohort->S_lamb[da];
                    for (a = 0; a <= sale_month; a++)
                        cohort->I[a] += cohort->I_lamb[da][a];
                }

                // sell oldest cohort
                cohorts[oldest_cohort].exists = false;
                oldest_cohort++;
            }
        }

        // Now at beginning of April, add a new cohort of lambs
        youngest_cohort++;
        lamb_cohort = &cohorts[youngest_cohort];
        set_cohort(lamb_cohort, 0, uninfected_cohort_size[0]);
        for (da = age_first_mating+1; da <= max_age; da++) {
            dam_cohort = get_cohort_of_age(da, cohorts, n_cohorts);
            dams = get_number_in_cohort_in_state(Any, dam_cohort, NULL);
            infectious_dams = get_number_in_cohort_in_state(Infectious, dam_cohort, p_infectious);
            lamb_cohort->I_lamb[da][0] = infectious_dams * lambs_per_dam[da] * prob_dam_to_lamb;
            lamb_cohort->S_lamb[da] = dams * lambs_per_dam[da] - lamb_cohort->I_lamb[da][0];
            lamb_cohort->S -= lamb_cohort->I_lamb[da][0];
            lamb_cohort->I[0] += lamb_cohort->I_lamb[da][0];
        }
    }

    // printf("\n ##### det simulation #####\n");
    // for (y = 0; y < sim_years; y++) {
    //     for (m = April; m <= March; m++) {
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
    //                 printf("% 7.4f ", r->infected_by_age[month][i]);
    //             printf("= ");
    //             printf("% 4.0f | ", r->infected[month]);
    //             for (i = 0; i <= max_age; i++)
    //                 printf("% 3.0f ", r->infectious_by_age[month][i]);
    //             printf("= ");
    //             printf("% 4.0f | ", r->infectious[month]);
    //             // for (i = 0; i <= max_age; i++)
    //             //     printf("% 5.0f ", r->total_by_age[month][i]);
    //             // printf("= ");
    //             // printf("% 6.0f | ", r->total[month]);
    //             printf("% 3.0f\n", r->flock_size[month]);
    //         }
    //     }
    // }

    cohorts_free(cohorts, n_cohorts);
    free(litter_dist_by_age[0]);
    free(litter_dist_by_age);
    free(Lambda);
    free(p_infectious);
    free(removal_rate);
    free(infection_rate);
    free(uninfected_cohort_size);
    free(contact_weights[0][0]);
    free(contact_weights[0]);
    free(contact_weights);
    return r;
}

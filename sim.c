d#include "R0.h"

double constant_inf(int a, params_t *p)
{
    return 1;
}

double linear_inf(int a, params_t *p)
{
    return a;
}

double exponential_inf(int a, params_t *p)
{
    return pow(10, (double)a/12.);
}

double four_yos_inf(int a, params_t *p)
{
    if (a < 48)
        return 0;
    else
        return 1;
}

void mix_in_field(int youngest, int oldest, double ***contact_weights, params_t *p) {
    int m, c, i;

    for (m = April; m <= March - p->housing_period; m++)
        for (i = youngest; i <= oldest; i++)
            for (c = i; c <= oldest; c++)
                contact_weights[m][i][c] = contact_weights[m][c][i] = 1;
}

void mix_in_housing(int youngest, int oldest, double pens, double ***contact_weights, params_t *p) {
    int m, c, i;

    for (m = March - p->housing_period + 1; m <= March; m++)
        for (i = youngest; i <= oldest; i++)
            for (c = i; c <= oldest; c++)
                contact_weights[m][i][c] = contact_weights[m][c][i] = 1/pens + p->inter_pen_contact_weight;
}

double*** all_mixed(params_t *p) {
    /*
        all ewes mix with all other ewes throughout the year + lambs (not reality)
    */
    int max_age = p->max_age;
    double pens = p->pens;
    double ***contact_weights = init_contact_weights(max_age);

    mix_in_field(0, max_age, contact_weights, p);
    mix_in_housing(0, max_age, pens, contact_weights, p);
    return contact_weights;
}

double*** mixed(params_t *p) {
    /*
        all ewes mix with all other ewes throughout the year
        lambs segregated from ewes and no contact between lambs 
        (perhaps not correct but prevents lambs infecting each other)
    */
    // int i, m;
    int max_age = p->max_age;
    double pens = p->pens;
    double ***contact_weights = init_contact_weights(max_age);

    // lambs have no contact with anything else
    // not needed as all weights are initialised to 0
    // for (m = April; m <= March; m++)
    //     for (i = 0; i <= max_age; i++)
    //             contact_weights[m][i][0] = contact_weights[m][0][i] = 0;

    mix_in_field(1, max_age, contact_weights, p); // ewes have field trans rate in field
    mix_in_housing(1, max_age, pens, contact_weights, p); // ewes have housing trans rate when housed
    return contact_weights;
}

double*** one_yos_seg_rest_mixed(params_t *p) {
    int max_age = p->max_age;
    double pens = p->pens;
    double ***contact_weights = init_contact_weights(max_age);

    mix_in_field(1, 1, contact_weights, p);
    mix_in_field(2, max_age, contact_weights, p);
    mix_in_housing(1, 1, 1, contact_weights, p);
    mix_in_housing(2, max_age-1, pens, contact_weights, p);
    return contact_weights;
}

double*** one_and_two_yos_mixed_rest_mixed(params_t *p) {
    int max_age = p->max_age;
    double pens = p->pens;
    double ***contact_weights = init_contact_weights(max_age);

    mix_in_field(1, max_age, contact_weights, p);
    mix_in_housing(1, 2, pens/2., contact_weights, p);
    mix_in_housing(3, max_age-1, pens/2., contact_weights, p);
    return contact_weights;
}

double*** oldest_seg_rest_mixed(params_t *p) {
    int max_age = p->max_age;
    double pens = p->pens;
    double ***contact_weights = init_contact_weights(max_age);

    mix_in_field(1, max_age-1, contact_weights, p);
    mix_in_housing(1, max_age-2, pens-1, contact_weights, p);
    mix_in_housing(max_age-1, max_age-1, 1, contact_weights, p);
    return contact_weights;
}

double*** youngest_oldest_seg_rest_mixed(params_t *p) {
    int max_age = p->max_age;
    double pens = p->pens;
    double ***contact_weights = init_contact_weights(max_age);

    mix_in_field(1, max_age-1, contact_weights, p);
    mix_in_housing(1, 1, 1, contact_weights, p);
    mix_in_housing(2, max_age-2, pens-2, contact_weights, p);
    mix_in_housing(max_age-1, max_age-1, 1, contact_weights, p);
    return contact_weights;
}

double*** all_segregated(params_t *p) {
    /*  
    Pasture
        all crops segregated
    Housing
        all crops segregated between pens with inter-pen transmission
    */
    int max_age = p->max_age;
    double ***contact_weights = init_contact_weights(max_age);
    return contact_weights;
}

void set_default_parameters(char *filename, char *text, int i, params_t *p) {
    strcpy(p->text, text);
    sprintf(p->file, "%s%d", filename, i);

    p->flock_size = 30;                // number of mated ewes in mating month
    p->pens = 1;                       // number of pens
    p->sim_years = 20;                 // number of years to simulate
    p->max_age = 5;                    // maximum age of a ewe in years
    p->age_first_mating = 1;           // age in years when first mated
    p->weaning_month = April;          // month that lambs are weaned
    p->sale_month = September;         // month that lambs and oldest ewes are sold
    p->mating_month = October;         // month in which the number of mated ewes equals flock_size
    p->housing_period = 3;             // number of months ewes are housed
    p->inter_pen_contact_weight = 0.;  // weighting of inter-pen transmission (0 to 1)
    p->sale_strategy = Random;         // how to sell lambs: Random or Oldest_dam_first

    p->L_mean = 15;                    // mean latent period
    p->L_sd = 5;                       // standard deviation in latent period
    p->S_mean = 8;                     // mean seroconversion period
    p->S_sd = 4;                       // standard deviation in seroconversion period
    // p->beta_field = 7.5e-4;         // per capita transmission rate when grazed
    p->beta_field = 0;                 // per capita transmission rate when grazed
    p->beta_housed = 0.17;             // per capita transmission rate when housed
    p->prob_dam_to_lamb = 0.1;         // probability an infectious dam infects each of its lambs
    p->prop_resistant = 0.0;           // proportion of ewes resistant
    p->tsi_removal = 1000;             // ewes with tsi greater than this value are sold
    p->gr_tsi_grad = 0;

    // beta_housed = 0.17 corresponds to 43 ewes housed under typical stocking density (~1 m2/ewe)
    // if we assume that stocking density is independent of flock size and contact rate linearly 
    // changes with flock size then increasing flock size means beta_housed must be changed to keep R0 constant

    // p->num_ibm_sims = 100;
    p->num_ibm_sims = 500;
    p->infectious_func = Constant;
    p->contact_weight_func = &mixed;
};

void MTvar3(params_t *p) {
    set_default_parameters("MTvar", "MTvar", 0, p);
    int i;
    int max_housed_days = 60;
    p->flock_size = 30;
    p->housing_period = max_housed_days / 30;
    p->num_ibm_sims = 200;
    p->sim_years = 300;
    double max_beta_housed = 0.17;
    double cdh[100];
    int n = 0;

    for (i = 0; i < 25; i++)
        cdh[n++] = i;
    for (i = 30; i <= max_housed_days; i+=10)
        cdh[n++] = i;

    p->prob_dam_to_lamb = 0;
    for (i = 0; i < n; i++) {
        p->beta_housed = max_beta_housed * cdh[i] / (double) max_housed_days;
        sprintf(p->file, "%s%d", "MTvar00", i);
        run_multiple(p);
    }

    p->prob_dam_to_lamb = 0.1;
    for (i = 0; i < n; i++) {
        p->beta_housed = max_beta_housed * cdh[i] / (double) max_housed_days;
        sprintf(p->file, "%s%d", "MTvar01", i);
        run_multiple(p);
    }

    p->prob_dam_to_lamb = 0.2;
    for (i = 0; i < n; i++) {
        p->beta_housed = max_beta_housed * cdh[i] / (double) max_housed_days;
        sprintf(p->file, "%s%d", "MTvar02", i);
        run_multiple(p);
    }
}

int main() {
    params_t p;
    set_seedRNG(0, &p);

    MTvar3(&p);
    
    gsl_rng_free(p.stream);
}


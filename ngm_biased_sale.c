#include "R0.h"

typedef struct {
    double eigenvalue;
    gsl_vector *eigenvector;
} eigen_t;

double prob_alive(int age, int age_at_infection, double *Lambda[], double *p_dam_age_before_sale, double *p_dam_age_after_sale, 
    int age_first_mating, int max_age, enum month sale_month) {
    int i;
    double p = 0;
    
    for (i = age_first_mating+1; i <= max_age; i++)
        if (isfinite(Lambda[i][age])) {
            if (age < sale_month)
                p += p_dam_age_before_sale[i] * exp(-Lambda[i][age] + Lambda[i][age_at_infection]);
            else
                p += p_dam_age_after_sale[i] * exp(-Lambda[i][age] + Lambda[i][age_at_infection]);
        }
    return p;
}

void dominant_eigen_solution(int n_age_classes, double M[][n_age_classes], eigen_t *eigen) {
    int j, k, status;
    gsl_matrix *A = gsl_matrix_alloc(n_age_classes, n_age_classes);
    gsl_vector_complex *eigenvalues = gsl_vector_complex_alloc(n_age_classes);
    gsl_matrix_complex *eigenvectors = gsl_matrix_complex_alloc(n_age_classes, n_age_classes);
    gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(n_age_classes);

    for (j = 0; j < n_age_classes; j++)
        for (k = 0; k < n_age_classes; k++)
            gsl_matrix_set(A, j, k, M[j][k]);

    gsl_set_error_handler_off();
    // gsl_eigen_nonsymm(A, eigenvalues, eigenvectors, w);
    // gsl_eigen_nonsymm_params(0, 1, w);
    status = gsl_eigen_nonsymm(A, eigenvalues, w);
    if (status) printf("Not all eigenvalues found\n");

    // view of real parts of vector of eigenvalues
    gsl_vector_view view_eigenvalues = gsl_vector_complex_real(eigenvalues);
    // get the maximum eigenvalue (R0)
    eigen->eigenvalue = gsl_vector_max(&view_eigenvalues.vector);

    // removed eigenvectors as not needed. to reinstate change nonsymm to nonsymmv in function calls
    // // get a view of the eigenvector associated with this eigenvalue
    // // index of maximum eigenvalue
    // size_t idx = gsl_vector_max_index(&view_eigenvalues.vector);
    // // a view of the complex eigenvector corresponding to the maximum eigenvalue
    // gsl_vector_complex_view view_eigenvector = gsl_matrix_complex_column(eigenvectors, idx);
    // // views of the real part of this complex eigenvector
    // gsl_vector_view view_eigenvector_real = gsl_vector_complex_real(&view_eigenvector.vector);
    // // copy the view of the real parts of the eigenvector into a new vector

    // gsl_vector_memcpy(eigen->eigenvector, &view_eigenvector_real.vector);
 
    gsl_eigen_nonsymm_free(w);
    gsl_matrix_complex_free(eigenvectors);
    gsl_vector_complex_free(eigenvalues);
    gsl_matrix_free(A);
}

void R0(int n_age_classes, params_t *p, double *lambs_per_dam, double *p_infectious, double *cohort_removal_rate, double *infection_rate, 
    double *uninfected_cohort_size, double ***contact_weights, double *Lambda, double *Lambda_age_0[], double *p_dam_age_before_sale, double *p_dam_age_after_sale, 
    eigen_t *eigen) {
    int max_age = p->max_age;
    int age_first_mating = p->age_first_mating;
    enum month weaning_month = p->weaning_month;
    enum month sale_month = p->sale_month;
    double prob_dam_to_lamb = p->prob_dam_to_lamb;
    double beta_housed = p->beta_housed;

    int c, i, a, m, tsi;
    int ae, ac;
    double foi, lamb_foi, p_alive, infectiousness;
    double M[n_age_classes][n_age_classes];

    // set matrix elements to zero
    for (a = 0; a < n_age_classes; a++)
        for (i = 0; i < n_age_classes; i++)
            M[a][i] = 0;

    // for a ewe infected at age a, find the expected number of ewes it 
    // infects of age i over its infected lifespan
    for (a = 0; a < n_age_classes; a++) {
        // for all cohorts from past to future that this infected ewe will come into contact with
        // +ve c are past cohorts
        // -ve c are future cohorts
        for (c = max_age; c >= -max_age; c--) {
            // ae: age of infected ewe when it first contacts the cohort
            // ac: age of cohort when it first contacts infected ewe
            if (a < -12*c) {
                // cohort is in the future so it is exposed from age 0
                // and the ewe is a multiple of 12 months old
                ae = -12*c;
                ac = 0;
            }
            else {
                // cohort is in the present or past so it is exposed 
                // from when this ewe becomes infected at age a
                // ae = age when this ewe can first infect susceptible ewes
                // ac = age of cohort when the ewe can first infect
                ae = a + 1;
                ac = ae + 12*c;
            }

            // number of susceptible ewes in this cohort the infected ewe infects over the cohorts lifespan
            while (ac < n_age_classes && ae < n_age_classes) {
                m = ae % 12;
                tsi = ae - a - 1;

                p_alive = prob_alive(ae-1, a, Lambda_age_0, p_dam_age_before_sale, p_dam_age_after_sale, 
                    age_first_mating, max_age, sale_month);

                // if lambs have just been born in this cohort and this ewe is a dam 
                // it contributes some infected lambs via infected colostrum
                if (prob_dam_to_lamb && ac == 0 && ae >= 12 * (age_first_mating+1) + April) {
                    M[a][0] += lambs_per_dam[ae/12] * prob_dam_to_lamb * p_infectious[tsi] * p_alive;
                    // if (a == 0)
                    //     printf("%f %f %f %f %f\n", M[a][0], lambs_per_dam[ae/12], prob_dam_to_lamb, p_infectious[tsi], p_alive);
                }
                infectiousness = p_infectious[tsi] * p->inf_fn_scale * (*p->inf_fn)(tsi, p);

                // force of infection exerted by this ewe
                foi = contact_weights[m][ae/12][ac/12] * infection_rate[m] * infectiousness;

                // maternal transmission via exhalation of virus during weaning
                // if (ac < weaning_month && ae >= 12 * (age_first_mating+1) + April)
                //     lamb_foi = lambs_per_dam[ae/12] * beta_housed * infectiousness;
                // else
                lamb_foi = 0;

                M[a][ac] += (uninfected_cohort_size[ac] * (1-exp(-foi)) + 1-exp(-lamb_foi)) * p_alive;

                ae++;
                ac++;
            }
        }
    }
    // for (a = 0; a < n_age_classes; a++)
    //     printf("%d\t%f\t%f\t%f\n", a, M[0][a]/M[1][a], M[0][a]/M[2][a], M[13][a]/M[14][a]);
    // double ss = 0;
    // for (a = 0; a < n_age_classes; a++) {
    //     printf("%d\t%f\t%f\n", a, M[a][0], M[a][16]);
    //     ss += M[a][a];
    // }
    // printf("ss = %f\n", ss);
    dominant_eigen_solution(n_age_classes, M, eigen);
    
    // for (a = 0; a < n_age_classes; a++) {
    //     printf("% 3d %.3lf  ", a, gsl_vector_get(eigen->eigenvector, a));
    //     for (i = 0; i < n_age_classes; i++)
    //         if (M[a][i] == 0)
    //             printf(".");
    //         else
    //             printf("*");
    //     printf("\n");
    // }
    // for (a = 0; a < n_age_classes; a++)
    //     printf("% 3d %.3lf\n", a, M[a][0]);

}

result_t* Next_Gen_Matrix_R0(params_t *p) {
    // calculate R0 from next generation matrix
    int max_age = p->max_age;
    int age_first_mating = p->age_first_mating;
    int n_age_classes = p->n_age_classes;
    enum month sale_month = p->sale_month;

    double *p_infectious = prob_infectious(p);
    double *base_removal_rate = natural_deaths(p->n_age_classes);
    double *infection_rate = monthly_infection_rate(p->housing_period, p->beta_housed, p->beta_field);
    double *uninfected_cohort_size = cohort_size(p, base_removal_rate);
    double ***contact_weights = (*p->contact_weight_func)(p);

    int i, a;
    double sum, *Lambda, *Lambda_age_0[max_age+1];
    double lambs_per_dam[max_age+1];
    double p_dam_age_before_sale[max_age+1];
    double p_dam_age_after_sale[max_age+1];
    double alpha[max_age+1];
    double ewe_removal_rate[n_age_classes];
    double cohort_removal_rate[n_age_classes];
    result_t *r = init_results(0, 0, n_age_classes);
    eigen_t eigen;
    eigen.eigenvector = gsl_vector_alloc(n_age_classes);

    double **litter_dist_by_age = litter_dist_by_age_data(age_first_mating, max_age);
    for (a = age_first_mating+1; a <= max_age; a++) {
        lambs_per_dam[a] = 0;
        for (i = 0; i < 5; i++)
            lambs_per_dam[a] += i * litter_dist_by_age[a][i];
    };
    // printf("---- lambs per dam by age\n");
    // for (a = age_first_mating+1; a <= max_age; a++)
    //     printf("%d %lf\n", a, lambs_per_dam[a]);

    // copy base removal rate into ewe and susceptible cohort removal rates
    for (i = 0; i < n_age_classes; i++) {
        cohort_removal_rate[i] = base_removal_rate[i];
        ewe_removal_rate[i] = base_removal_rate[i];
    }

    // calculate age-specific retain probabilities TODO MOVE THIS OUT ALSO IN DET
    get_retain_probabilities(alpha, uninfected_cohort_size, lambs_per_dam, p);
    // printf("---- retain probabilities\n");
    // for (i = age_first_mating+1; i <= max_age; i++)
    //     printf("%d %lf\n", i, alpha[i]);

    // calculate removal rate of sold lambs in susceptible cohort
    cohort_removal_rate[sale_month-1] = -log(uninfected_cohort_size[sale_month] / uninfected_cohort_size[sale_month-1]);
    Lambda = cumulative_removal_rate(n_age_classes, cohort_removal_rate);
    // printf("---- cohort removal rate\n");
    // for (i = 0; i < n_age_classes; i++)
    //     printf("%d %lf\n", i, cohort_removal_rate[i]);


    // total number of dams at beginning of April
    sum = size(April, age_first_mating+1, max_age, n_age_classes, uninfected_cohort_size);
    // printf("dams in April = %lf\n", sum);

    for (i = age_first_mating+1; i <= max_age; i++) {
        // calculate probability that the infected ewe has a dam of age a
        p_dam_age_before_sale[i] = uninfected_cohort_size[12*i + April] / sum;
        // adjust ewe removal rate on sale month
        ewe_removal_rate[sale_month-1] = -log(alpha[i]);
        // survival vectors for ewes of each dam age
        Lambda_age_0[i] = cumulative_removal_rate(n_age_classes, ewe_removal_rate);
    }
    // printf("---- ewe removal rate\n");
    // for (i = 0; i < n_age_classes; i++)
    //     printf("%d %f\n", i, ewe_removal_rate[i]);

    sum = 0;
    for (i = age_first_mating+1; i <= max_age; i++) {
        p_dam_age_after_sale[i] = alpha[i] * p_dam_age_before_sale[i];
        sum += p_dam_age_after_sale[i];
    }
    for (i = age_first_mating+1; i <= max_age; i++)
        p_dam_age_after_sale[i] /= sum;
    // for (i = age_first_mating+1; i <= max_age; i++)
    //     printf("%d %lf %lf\n", i, p_dam_age_before_sale[i], p_dam_age_after_sale[i]);

    // calculate R0
    R0(n_age_classes, p, lambs_per_dam, p_infectious, cohort_removal_rate, infection_rate, 
        uninfected_cohort_size, contact_weights, Lambda, Lambda_age_0, p_dam_age_before_sale, p_dam_age_after_sale, &eigen);

    // save relevant values for return in r
    r->R0 = eigen.eigenvalue;
    r->generation_time = mean_infected_life_expectancy(p, uninfected_cohort_size, Lambda);
    r->doubling_time = log(2)/log(r->R0) * r->generation_time;
    r->init_time = p->mating_month;

    // free memory
    gsl_vector_free(eigen.eigenvector);
    free(litter_dist_by_age[0]);
    free(litter_dist_by_age);
    free(Lambda);
    for (i = age_first_mating+1; i <= max_age; i++)
        free(Lambda_age_0[i]);
    free(p_infectious);
    free(base_removal_rate);
    free(infection_rate);
    free(uninfected_cohort_size);
    free(contact_weights[0][0]);
    free(contact_weights[0]);
    free(contact_weights);
    return r;
}


import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from cycler import cycler

color_cycle = cycler('color', ['C0', 'C1', 'C2'])
marker_cycle = cycler('marker', ['o', 's', '^'])
linestyle_cycle = cycler('ls', ['-', '-', '-'])
hatch_cycle = cycler('hatch', ['-', '\\', '/'])

def new_dict(n_cohorts):
    x = [[] for _ in range(n_cohorts)]
    return {'year':[],
        'total_infected':[],
        'total_positive':[],
        'total_resistant':[],
        'seroconversions':[],
        'prop_pos':[],
        'num_infected_in_flock':[],
        'num_positive_in_flock':[],
        'flock_infected_prevalence':[],
        'flock_sero_prevalence':[],
        'mean_infectious_age':[],
        'sero_prevalence':deepcopy(x),
        'prevalence':deepcopy(x),
        'infectiousness':deepcopy(x),
        'foi':deepcopy(x),
        'total':deepcopy(x),
        'positive':deepcopy(x),
        'infected':deepcopy(x),
        'tsi':deepcopy(x),
    }

def get_results(filename):
    ngm = {}
    parameters = {}
    stochastic = []
    deterministic = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith("stochastic"):
                n_cohorts = int(line.split()[-1])
                stochastic.append(new_dict(n_cohorts))
                data = stochastic[-1]

            elif line.startswith("deterministic"):
                n_cohorts = int(line.split()[-1])
                deterministic = new_dict(n_cohorts)
                data = deterministic

            elif line.startswith("matrix"):
                ngm['R0'] = float(next(f).split()[1])
                ngm['gen_time'] = float(next(f).split()[1])/12
                ngm['double_time'] = float(next(f).split()[1])/12
                ngm['init_time'] = float(next(f).split()[1])/12

            elif line.startswith("parameters"):
                parameters['sim_years'] = int(next(f).split()[1])
                parameters['flock_size'] = float(next(f).split()[1])
                parameters['all_ewes_flock_size'] = float(next(f).split()[1])
                parameters['max_age'] = int(next(f).split()[1])
                parameters['age_first_mating'] = int(next(f).split()[1])
                parameters['pens'] = float(next(f).split()[1])
                parameters['housing_period'] = int(next(f).split()[1])
                parameters['beta_housed'] = float(next(f).split()[1])
                parameters['beta_field'] = float(next(f).split()[1])
                parameters['L_mean'] = round(float(next(f).split()[1]))
                parameters['L_sd'] = round(float(next(f).split()[1]))
                parameters['S_mean'] = round(float(next(f).split()[1]))
                parameters['S_sd'] = round(float(next(f).split()[1]))
                parameters['inter_pen_contact_weight'] = float(next(f).split()[1])
                parameters['prob_dam_to_lamb'] = float(next(f).split()[1])
                parameters['text'] = next(f).rstrip()

            else:
                x = iter([float(i) for i in line.split()])
                for y in x:
                    if y % 12 != 11:
                        break
                    data['year'].append(y/12)
                    data['flock_infected_prevalence'].append(next(x))
                    
    return parameters, stochastic, deterministic, ngm

def mean_infected(filenames):
    fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8))
    labels = ['0.0', '0.1', '0.2']
    pmt = np.array(list(range(25)) + list(range(30, 70, 10)))

    # R0 taken from output of simulations for MT=0
    R0 = [
        0.122642, 0.245119, 0.367433, 0.489583,
        0.611569, 0.733392, 0.855052, 0.976549,
        1.097884, 1.219056, 1.340067, 1.460915,
        1.581602, 1.702128, 1.822493, 1.942697,
        2.062740, 2.182623, 2.302346, 2.421909,
        2.541312, 2.660557, 2.779642, 2.898568,
        3.017335, 3.726623, 4.896240, 6.050398,
        7.189333,
    ]
    n = len(pmt)

    for filename, label, m in zip(filenames[::-1], labels[::-1], (color_cycle+hatch_cycle)[::-1]):
        msp = []
        lsp = []
        usp = []
        for i in range(n):
            f = filename + str(i) + '.txt'
            results = get_results(f'Results/{f}')
            # all stochastic sims at this R0
            stochastic = results[1]

            # array of mean prevalences, one for each simulation 
            d = 100*np.array([np.mean(z['flock_infected_prevalence'][-100:]) for z in stochastic])
            # remove sims that never took off
            d = d[d>0]
            # but if there are no sims that took off then add one zero to d
            if len(d) == 0:
                d = np.zeros(1)
            # the mean equilibrium prevalence over all sims
            mean_equilibrium_prevalence = np.mean(d)
            msp.append(mean_equilibrium_prevalence)
        
        ax.plot(R0, msp, color=m['color'], lw=0.5, label=label)

    ax.set_ylabel('Mean 10-year prevalence\n(percentage of ewes infected)')
    ax.axvline(1, ls=':', color='k')
    ax.set_xlabel('R0 (horizontal transmission only)')
    plt.legend(title='Maternal transmission rate', loc='lower right')

    fig.savefig('fig1.png')

mean_infected(['MTvar00', 'MTvar01', 'MTvar02', ])

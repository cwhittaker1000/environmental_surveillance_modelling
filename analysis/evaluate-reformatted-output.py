# Usage:
#    $ Rscript analysis/dump-arriving-reads.R | \
#        python3 analysis/reformat-r-output.py | \
#        python3 analysis/evaluate-reformatted-output.py
#    ...
#    0: detected 9.8% daily growth at p=3.55e-11 on day=114, when 49% had been
#    infected in A

import sys
import math
import statsmodels.api as sm

MIN_WINDOW_SIZE=60
POPULATION_A=10**7
INITIAL_INFECTIONS=10

days = []
read_simulations = []
S_simulations = []
for line in sys.stdin:
    line = line.strip()
    day, *vals = line.split('\t')
    days.append(int(day))

    if not read_simulations:
        for i in range(len(vals)):
            read_simulations.append([])
            S_simulations.append([])

    for n_simulation, reads_S in enumerate(vals):
        reads, S = reads_S.split(":")
        read_simulations[n_simulation].append(int(reads))
        S_simulations[n_simulation].append(int(S))

if S_simulations[0][0] != POPULATION_A - INITIAL_INFECTIONS:
    raise Exception("Initial conditions out of sync with input")


for n_simulation, read_simulation in enumerate(read_simulations):
    best_day = None
    best_pvalue = None
    best_growth = None
    for day in range(MIN_WINDOW_SIZE, max(days)):
        pvalue = 1
        coef = 0

        read_simulation_subset = read_simulation[:day]
        days_subset = days[:day]

        # can't possibly detect exponential growth, and p-values are not
        # reliable for things like 0 0 0 ... 1
        if sum(read_simulation_subset) <= 5: continue

        model = sm.GLM(read_simulation_subset, sm.add_constant(days_subset),
                       family=sm.families.Poisson())
        try:
            result = model.fit()
        except ValueError:
            continue

        pvalue = result.pvalues[1]
        coef = result.params[1]

        # In log space; convert to percentage daily growth.
        coef = math.exp(coef)-1

        # Skip things that are decreasing
        if coef < 0:
            continue

        if best_day is None or pvalue < best_pvalue:
            best_day = day
            best_pvalue = pvalue
            best_growth = coef

            if best_pvalue < 1e-10:
                #print(" ".join(str(x) for x in days_subset))
                #print(" ".join(str(x) for x in read_simulation_subset))

                S = S_simulations[n_simulation][best_day]
                cumulatively_infected = POPULATION_A - S

                print("%s: detected %.1f%% daily growth at p=%.2e on day=%s, "
                      "when %.0f%% had been infected in A" % (
                          n_simulation, best_growth*100, best_pvalue, best_day,
                          cumulatively_infected *  100 / POPULATION_A
                      ))
                break

# Usage:
#    $ Rscript analysis/dump-arriving-reads.R | \
#        python3 analysis/reformat-r-output.py | \
#        python3 analysis/evaluate-reformatted-output.py
#    ...
#    0: detected 2.1% daily growth at p=8.37e-11 on day=316

import sys
import math
import statsmodels.api as sm

MIN_WINDOW_SIZE=60

days = []
simulations = []
for line in sys.stdin:
    line = line.strip()
    day, *vals = line.split('\t')
    days.append(int(day))

    if not simulations:
        for i in range(len(vals)):
            simulations.append([])

    for n_simulation, val in enumerate(vals):
        simulations[n_simulation].append(int(val))
    

for n_simulation, simulation in enumerate(simulations):
    best_day = None
    best_pvalue = None
    best_growth = None
    for day in range(MIN_WINDOW_SIZE, max(days)):
        pvalue = 1
        coef = 0

        simulation_subset = simulation[:day]
        days_subset = days[:day]

        # can't possibly detect exponential growth, and p-values are not
        # reliable for things like 0 0 0 ... 1
        if sum(simulation_subset) <= 5: continue
        
        model = sm.GLM(simulation_subset, sm.add_constant(days_subset),
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
                #print(" ".join(str(x) for x in simulation_subset))
                print("%s: detected %.1f%% daily growth at p=%.2e on day=%s" % (
                    n_simulation, best_growth*100, best_pvalue, best_day))
                break



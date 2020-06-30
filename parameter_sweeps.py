#!/usr/bin/env python

# Script to setup parameter sweeps over different sets of parameters
# and get corresponding plots

import os, sys
from random import shuffle
import numpy as np
from pickle import load, dump
import random
from combined_simulations import Get_Abundance_Data, core_sim_loop, average_simulation_data
from protein_abundances import generate_complex_ids, last_complex_adjustment, align_complex_abundances, sorted_complex_abundances


alpha_factor = 1.0
alpha=1
ideality_correction=1
deviation_from_ideality = 1.0
abundance_correlation = 0.7
large_complex_boost = 20
large_complex_correlation = 0.85

Sweep_Parameters = {
    "complex size" : {
        "sweep_elem" : "total_partners",
        "total_partners" : [2, 5, 10, 20, 40],
        "abundance_correlation" : 0.7,
        "alpha" : 1,
        "ideality_correction" : 1
        },
    "abundance correlation" : {
        "sweep_elem" : "abundance_correlation",
        "total_partners" : 20,
        "abundance_correlation" : np.linspace(0.5, 0.9, 5).tolist(),
        "alpha" : 1,
        "ideality_correction" : 1
        },
    "water abundance" : {
        "sweep_elem" : "alpha",
        "total_partners" : 20,
        "abundance_correlation" : 0.7,
        "alpha" : np.linspace(0.15, 1.0, 6).tolist(),
        "ideality_correction" : 1
        },
    "ideality correction" : {
        "sweep_elem" : "ideality_correction",
        "total_partners" : 20,
        "abundance_correlation" : 0.7,
        "alpha" : 1,
        "ideality_correction" : np.linspace(0.5, 1.5, 6).tolist()
        }
    }

def sweep_parameter(base, sweep_name):
    """
    Performs a sweep over complex sizes
    """
    sweep_key = Sweep_Parameters[sweep_name]["sweep_elem"]
    sweep_values = Sweep_Parameters[sweep_name][sweep_key]
    Data = { key : [] for key in sweep_values }
    arr_base = np.array(base) + 1
    Sweepfile = "%s.dmp" % sweep_name.capitalize().replace(" ","_")

    for key in sweep_values:

        Sweep_Parameters[sweep_name][sweep_key] = key
        complex_contents = generate_complex_ids([Sweep_Parameters[sweep_name]["total_partners"]], len(abundance_range))
        complex_contents = last_complex_adjustment(complex_contents, len(abundance_range))
        aligned_abundances = align_complex_abundances(complex_contents, abundance_range, abundance_correlation = Sweep_Parameters[sweep_name]["abundance_correlation"])
        aligned_abundances = sorted_complex_abundances(aligned_abundances, complex_contents[-1], abundance_range, abundance_correlation = Sweep_Parameters[sweep_name]["abundance_correlation"])
        re_runs, buckets = core_sim_loop(base, complex_contents, aligned_abundances)
        means, stds, pre_buckets = average_simulation_data(re_runs, buckets, alpha = Sweep_Parameters[sweep_name]["alpha"], ideality_correction = Sweep_Parameters[sweep_name]["ideality_correction"])
        Data[key] = (means, stds, pre_buckets)

    with open(Sweepfile, "w") as fp: dump({"Sweep_Data" : Data, "arr_base" : arr_base}, fp)

if __name__ == '__main__':
    
    Datastatfile = "data_stats_dump.dmp"
    Interaction = "Paxdb"

    abundance_range, total_partners = Get_Abundance_Data(Datastatfile, Interaction = Interaction)
    base = np.linspace(0.0, 1.0, 20).tolist()
    arr_base = np.array(base) + 1
    sweep_parameter(base, "complex size")
    sweep_parameter(base, "abundance correlation")
    sweep_parameter(base, "water abundance")
    sweep_parameter(base, "ideality correction")

    #total_partners = [20]
    #abundance_correlation_list = np.linspace(0.5, 0.9, 5).tolist()
    #sweep_abundance_correlation(base, abundance_correlation_list)

    #total_partners = [20]
    #alpha_list = np.linspace(0.15, 0.95, 5).tolist()
    #sweep_abundance_correlation(base, alpha_list)

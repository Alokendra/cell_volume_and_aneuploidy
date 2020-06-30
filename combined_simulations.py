#!/usr/bin/env python

import os, sys
from random import shuffle
import numpy as np
from pickle import load, dump
import random
from protein_abundance_preprocess import Partners_Paxdb, Partners_Gamma
from protein_abundances import generate_complex_ids, last_complex_adjustment, align_complex_abundances, sorted_complex_abundances, calculate_free_mol_entities

random.seed(42)

alpha_factor = 1.0
alpha=1
ideality_correction=1
deviation_from_ideality = 1.0
abundance_correlation = 0.7
large_complex_boost = 20
large_complex_correlation = 0.85

def Get_Abundance_Data(Datastatfile, Interaction = "Paxdb"):
    """
    Obtains the interaction data from different database
    files
    """

    if Interaction == "Paxdb":
        abundance_range, total_partners = Partners_Paxdb(Datastatfile)
    elif Interaction == "Gamma":

        with open(Datastatfile) as fp:
            abundance_range, total_partners = load(fp)
        Gamma_Parameters = {
            "shape" : 1.538,
            "scale" : 2.016,
            "size" : len(total_partners)
        }
        abundance_range, total_partners = Partners_Gamma(Gamma_Parameters)
    else:
        raise Exception("Unrecognized interaction %s" % Interaction)

    return abundance_range, total_partners

def base_simulation(base, complex_contents, aligned_abundances, buckets):
    """
    Runs the simulation for a single replicate
    over the list of partners
    """
    read_out = []
    sub_bucket_list = []
    for aneuploidy_factor in base:
        sub_buckets = {}

        if buckets:
            for key in buckets.keys():
                sub_buckets[key] = [0, 0, 0]
        total_molecules, sub_buckets = calculate_free_mol_entities(aneuploidy_factor, complex_contents, aligned_abundances, sub_buckets)
        read_out.append(total_molecules)
        sub_bucket_list.append(sub_buckets)

    return read_out, sub_bucket_list


def core_sim_loop(base, complex_contents, aligned_abundances, repeats=5, buckets=[]):

    re_runs = []
    buckets = { key : [[]for _ in range(repeats)] for key in buckets }

    for i in range(0, repeats):
        read_out, sub_bucket_list = base_simulation(base, complex_contents, aligned_abundances, buckets)
        re_runs.append(np.array(read_out))
        for key in buckets.keys():
            buckets[key][i].append([elem[key] for elem in sub_bucket_list])

    re_runs = np.array(re_runs)

    return re_runs, buckets

def average_simulation_data(re_runs, buckets, alpha=1, ideality_correction=1):

    means = np.mean(re_runs, 0)
    Norm_Factor = (alpha*ideality_correction)/means[0]
    means = means * Norm_Factor
    stds = np.std(re_runs, 0) * Norm_Factor

    for key, value in buckets.iteritems():
        val = np.array(value)
        m = np.mean(val, axis = 0)
        m = m.reshape(-1, m.shape[-1])
        s = np.std(val, axis = 0)
        s = s.reshape(-1, s.shape[-1])
        buckets[key] = np.hstack((m, s))

    return means, stds, buckets


if __name__ == "__main__":

    Datastatfile = "data_stats_dump.dmp"
    Interaction = "Paxdb"

    abundance_range, total_partners = Get_Abundance_Data(Datastatfile, Interaction = Interaction)
    base = np.linspace(0.0, 1.0, 20).tolist()
    arr_base = np.array(base) + 1
    #sys.exit(1)
    complex_contents = generate_complex_ids(total_partners, len(abundance_range))
    complex_contents = last_complex_adjustment(complex_contents, len(abundance_range))
    aligned_abundances = align_complex_abundances(complex_contents, abundance_range, abundance_correlation = abundance_correlation)
    aligned_abundances = sorted_complex_abundances(aligned_abundances, complex_contents[-1], abundance_range, abundance_correlation =  abundance_correlation)
    re_runs, buckets = core_sim_loop(base, complex_contents, aligned_abundances, buckets = [3, 15])
    means, stds, pre_buckets = average_simulation_data(re_runs, buckets, alpha = alpha, ideality_correction = ideality_correction)

    Sim_Data = {"arr_base" : arr_base, "means" : means, "stds" : stds, "buckets" : pre_buckets}

    with open("simulation_data.dmp", "w") as fp: dump(Sim_Data, fp)


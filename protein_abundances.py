#!/usr/bin/env python

import os
import random
from random import shuffle
import numpy as np
from pickle import load


def calculate_van_Hoeff(aneuploidy_factor, complex_size):
    # We assume we have 1000 complexes that absorb 6000 genes total
    genes = range(0, 6000)
    # complexes are consecutive genes runs (simplifying assumption)
    shuffle(genes)

    genes = genes[:int(len(genes)*aneuploidy_factor)]
    genes = sorted(genes)
    genes = np.array(genes)

    complex_counter = 0
    free_protein_counter = len(genes)
    for a, b in zip(range(0, 6000-complex_size, complex_size), range(complex_size, 6000, complex_size)):
        if np.sum(np.logical_and(genes > a, genes < b))+1 == complex_size:
            complex_counter += 1
            free_protein_counter -= complex_size

    total_van_Hoeff = complex_counter + free_protein_counter

    return total_van_Hoeff

def generate_complex_ids(total_partners, abundance_range_length, manual_complex_size = 45):

    complex_contents = []
    i = 0
    local_partners = (np.array(total_partners)+1).tolist()

    while i < abundance_range_length - manual_complex_size:
        current_complex_size = random.choice(local_partners)
        complexes = []
        for j in range(i, i+current_complex_size):
            complexes.append(j)
        complex_contents.append(complexes)
        # print complex
        i += current_complex_size

    return complex_contents

def last_complex_adjustment(complex_contents, abundance_range_length):

    last_complex = complex_contents[-1]
    last_complex = np.array(last_complex[last_complex < abundance_range_length]).tolist()
    if type(last_complex) is int:
        del complex_contents[-1]
    else:
        complex_contents[-1] = last_complex

    # Manual large, abundant complex injection:
    complex_contents.append(range(complex_contents[-1][-1], abundance_range_length))

    return complex_contents


def align_complex_abundances(complex_contents, abundance_range, abundance_correlation=0.7, large_complex_boost = 20, large_complex_correlation = 0.85):

    aligned_abundances = np.copy(abundance_range)

    for complexes in complex_contents[:-1]:

        # Manual injection boosting of small complexes abundances:
        if len(complexes) < 45 and len(complexes) > 40:
            aligned_abundances[complexes] *= large_complex_boost
            loc_ab_corr = large_complex_correlation
        else:
            loc_ab_corr = abundance_correlation

        average_abundance = np.mean(aligned_abundances[complexes])
        aligned_abundances[complexes] = aligned_abundances[complexes]*(1 - loc_ab_corr) + average_abundance*loc_ab_corr

    return aligned_abundances

def sorted_complex_abundances(aligned_abundances, last_complex, abundance_range, abundance_correlation = 0.7):

    # Manual large, abundant complex injection:
    sorted_abundances = np.sort(abundance_range)
    average_abundance = np.mean(sorted_abundances[last_complex])
    aligned_abundances[last_complex] = sorted_abundances[last_complex]*(1-abundance_correlation) + average_abundance*abundance_correlation

    return aligned_abundances


def calculate_free_mol_entities(aneuploidy_factor, complex_contents, abundances, buckets={}):
    multiplication_factor = np.random.binomial(1, aneuploidy_factor, len(abundances))+1
    aneup_abundance = abundances * multiplication_factor

    total_proteins = sum(aneup_abundance)
    total_complexes = 0
    bound_proteins = 0

    for complexes in complex_contents:
        complex_abundances = aneup_abundance[complexes]
        loc_min = np.min(complex_abundances)
        total_complexes += loc_min
        bound_proteins += loc_min * len(complexes)

        complex_size = len(complexes)
        # yup, we are doing it with pointers and pointers are uncool and unpythonic.
        if buckets and complex_size in buckets.keys():
            buckets[complex_size][0] += loc_min
            buckets[complex_size][1] += np.sum(aneup_abundance[complexes])
            buckets[complex_size][2] += np.sum(aneup_abundance[complexes] - loc_min)

    total_molecules = total_proteins - bound_proteins + total_complexes

    return total_molecules, buckets

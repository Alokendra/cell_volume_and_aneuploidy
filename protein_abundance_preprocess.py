#!/usr/bin/env python

from pickle import load
import numpy as np

def Partners_Paxdb(Datastatfile, lower = 1, upper = 45):
    """
    Obtain the total number of partners from the
    PaxDB file
    """
    with open(Datastatfile) as fp:
        abundance_range, total_partners = load(fp)

    abundance_range = abundance_range[abundance_range > 1]
    total_partners = np.array(total_partners)
    total_partners = total_partners[total_partners > lower]
    total_partners = total_partners[total_partners < upper]
    total_partners = total_partners.tolist()

    return abundance_range, total_partners

def Partners_Gamma(Gamma_Parameters, lower = 1, upper = 45):
    """
    Obtain the number of partners from gamma distribution
    """
    
    total_partners = np.random.gamma(**Gamma_Parameters)
    total_partners = np.ceil(total_partners).astype(np.int16)
    total_partners = total_partners[total_partners > lower]
    total_partners = total_partners[total_partners < upper]
    total_partners_new = total_partners
    
    abundance_range = np.random.permutation(abundance_range)

    return abundance_range, total_partners

def Ploidy_Data(Ploidyfile):
    """
    Reads and transforms the ploidy data
    """
    with open(Ploidyfile) as fp: ploidy_vs_size = load(fp)

    corrfactor = ploidy_vs_size[0, 1]
    ploidy_vs_size[:, 1] /= corrfactor
    ploidy_vs_size[:, 2] /= corrfactor

    return ploidy_vs_size, corrfactor


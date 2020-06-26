#!/usr/bin/env python

"""
This script reads data from the interaction database and
ploidy vs size data, plots the density function and dumps
them as numpy array in separate files.
"""

from csv import reader as csv_reader
import numpy as np
from matplotlib import pyplot as plt
from collections import defaultdict
import os, pickle
from scipy.stats import gaussian_kde
from scipy import stats


def create_abundance_table(filename):
    """
    Creates the table of abundances from
    the input text file
    """
    abundance_range = []
    abundance_table = defaultdict(int)

    with open(filename) as source:
        reader = csv_reader(source, delimiter='\t')
        for line in reader:
            abundance = float(line[2])
            abundance_range.append(abundance)
            name = line[1].split('.')[1]
            abundance_table[name] = abundance

    return abundance_range, abundance_table


def create_interaction_table(filename):
    """
    Creates the interaction table from the interaction
    table database
    
    Parameters
    ----------
    `filename` : `str` Name of the input file
    
    Returns
    -------
    `list`
    A list containing names of interacting proteins
    """
    interaction_table = []
    with open(filename) as source:
        reader = csv_reader(source, delimiter='\t')
        header = reader.next()
        for line in reader:
            prot1 = line[0]
            prot2 = line[1]
            interaction_table.append([prot1, prot2])

    return interaction_table


def get_ploidy_data(filename):
    """
    Gets the ploidy data from the ploidy size file
    """
    exp_ploid_vs_size = []
    with open(filename, 'rb') as source:
        reader = csv_reader(source)
        reader.next()
        for line in reader:
            # print line
            ploidy = float(line[1])
            abs_radius = float(line[2])
            rel_radius = float(line[3])
            # print ploidy, rel_radius
            exp_ploid_vs_size.append([ploidy, abs_radius, rel_radius])

    return exp_ploid_vs_size


def Partner_Abundance(interaction_table):
    """
    Creates the list of partner abundances
    """
    

    total_partners = []
    partner_abundances = []
    
    for unique_gene in np.unique(interaction_table):
        c1 = np.count_nonzero(interaction_table[:, 0] == unique_gene)
        c2 = np.count_nonzero(interaction_table[:, 1] == unique_gene)
        interacting_partners = set(interaction_table[interaction_table[:, 0] == unique_gene, 1]) | set(interaction_table[interaction_table[:, 1] == unique_gene, 0])
        interacting_partners = list(interacting_partners)
        total_partners.append(len(interacting_partners))
        for partner in interacting_partners:
            partner_abundances.append([abundance_table[unique_gene], abundance_table[partner]])

    return total_partners, partner_abundances


def plot_hist():
    """
    Plots the histogram for the ranges
    """

    plt.hist(log_range, bins=50)
    plt.show()
    
    plt.hist(np.log(total_partners), bins=50)
    plt.show()
    
    plt.loglog(partner_abundances[:, 0], partner_abundances[:, 1], '.k')
    plt.show()


def plot_abundance_distribution(data):
    """
    Plots the abundance distribution of proteins
    """

    plt.title('Protein abundance distribution (log)')
    density = gaussian_kde(data.flatten())
    xs = np.linspace(data.min(), data.max(), 50)
    plt.plot(xs, density(xs), 'k')
    
    plt.xlabel('complex abundance (log, arbitrary units)')
    plt.ylabel('distribution density')
    plt.legend()
    plt.savefig("figures/Complex_Abundances.png")
    plt.close()

if __name__ == '__main__':
    
    Abundance_File = os.path.join("data", "4932-WHOLE_ORGANISM-integrated.txt")
    Interaction_File = os.path.join("data", "CervBinaryHQ.txt")
    Ploidy_File = os.path.join("data", "ploidy_size.csv")

    abundance_range, abundance_table = create_abundance_table(Abundance_File)
    abundance_range = np.array(abundance_range)
    log_range = np.log(abundance_range[abundance_range > 0])

    interaction_table = create_interaction_table(Interaction_File)
    interaction_table = np.array(interaction_table)

    exp_ploid_vs_size = get_ploidy_data(Ploidy_File)
    exp_ploid_vs_size = np.array(exp_ploid_vs_size)

    total_partners, partner_abundances = Partner_Abundance(interaction_table)
    partner_abundances = np.array(partner_abundances)

    #plot_hist()
    plot_abundance_distribution(log_range)

    slope, intercept, r_value, p_value, std_err = stats.linregress(partner_abundances[:, 0], partner_abundances[:, 1])
    print slope, intercept, r_value, p_value, std_err
    
    sel = np.logical_and(partner_abundances[:, 0] > 0, partner_abundances[:, 1] > 0)
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(partner_abundances[sel, 0]), np.log(partner_abundances[sel, 1]))
    print slope, intercept, r_value, p_value, std_err
    
    with open('ploidy_vs_size.dmp', 'w') as fp: pickle.dump(exp_ploid_vs_size, fp)
    with open('data_stats_dump.dmp', 'w') as fp: pickle.dump((abundance_range, total_partners), fp)

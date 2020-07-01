#!/usr/bin/env python

import os
from random import shuffle
import numpy as np
from matplotlib import pyplot as plt
from pickle import load
import random
from matplotlib.cm import get_cmap
#from statsmodels.nonparametric.smoothers_lowess import lowess
from calc_vol_press import Cell_Volumes, Osmotic_Pressure
from protein_abundance_preprocess import Ploidy_Data

Plotdir = "figures"

Datadir = "data"
Ploidyfile = os.path.join(Datadir, "ploidy_vs_size.dmp")
ploidy_vs_size, corrfactor = Ploidy_Data(Ploidyfile)
abundance_correlation = 0.7
large_complex_boost = 20
large_complex_correlation = 0.85

Suptitle1 = "large complex boost: x%.2f; large complex corr: %.2f, overall corr: %.2f\n" % (large_complex_boost, large_complex_correlation, abundance_correlation)

Base_Labels = { "complex size" : 1,
                "abundance correlation" : 0,
                "water abundance" : "base",
                "ideality correction" : "base" }

Label_Factor = { "complex size" : 1,
                 "abundance correlation" : 100,
                 "water abundance" : 100,
                 "ideality correction" : 100 }

def complex_size_swipe(Sweep_Data, arr_base, Plot_Parameters):
    """
    Plots a distribution of complex sizes
    """
    Title = Plot_Parameters["Title"]
    cmap = get_cmap("brg", len(Sweep_Data))
    plt.title(Title)
    plt.plot(arr_base, corrfactor*np.cbrt(arr_base), '--k', label=str(Plot_Parameters["base_label"]), lw=2)
    
    ax = plt.gca()
    ax.set_axisbelow(True)
    
    ax.yaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    ax.xaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    
    #plt.xticks(np.linspace(1, 2, 9), ['1.00', '', '1.25', '', '1.50', '', '1.75', '', '2.00'])
    #plt.yticks(np.linspace(8, 15, 8), [8, '', 10, '', 12, '', 14, ''])
    
    plt.xlabel("Ploidy")
    plt.ylabel("Equivalent Diameter")
    
    for ii, c_size in enumerate(sorted(Sweep_Data.keys())):
        c = cmap(ii)
    
        means, stds, buckets = Sweep_Data[c_size]
    
        plt.plot(arr_base,
                 corrfactor*np.cbrt(means),
                 '--',
                 color=c,
                 label=c_size * Plot_Parameters["label_factor"], lw=2)
    
        plt.fill_between(arr_base,
                         corrfactor*np.cbrt(means-stds),
                         corrfactor*np.cbrt(means+stds),
                         color=c,
                         alpha=.3)
    
    
    #plt.legend(loc='best', ncol=3, title='complex size')
    plt.legend(loc='best', ncol=3, title=Plot_Parameters["legend_title"])
    #plt.axis([0.95, 2.05, 8, 15])
    #plt.savefig(os.path.join(Plotdir, "%s.png" % Title.replace(" ", "_")))
    plt.savefig(Plot_Parameters["Figure_File_Name"])
    plt.close()

def plot_abundances(buckets):
    """
    Plots the abundance values of the proteins
    """
    for key, value_matrix in buckets.iteritems():
    
        plt.figure()
        plt.title('Molecules abundance vs ploidy for complex size %s' % key)
    
        #print value_matrix
    
        norm = value_matrix[0, 0] + value_matrix[0, 2]
        value_matrix /= norm
    
        plt.plot(arr_base, value_matrix[:, 0], '-k', label='complexes', lw=2)
    
        plt.plot(arr_base, value_matrix[:, 2], '-b', label='free proteins', lw=2)
    
        plt.plot(arr_base, value_matrix[:, 0] + value_matrix[:, 2], '-r', label='total molecules', lw=2)
    
        plt.fill_between(arr_base,
                         value_matrix[:, 0]+value_matrix[:, 3],
                         value_matrix[:, 0]-value_matrix[:, 3],
                         color='k',
                         alpha=.3)
    
        plt.fill_between(arr_base,
                         value_matrix[:, 2]+value_matrix[:, 5],
                         value_matrix[:, 2]-value_matrix[:, 5],
                         color='b',
                         alpha=.3)
    
        plt.fill_between(arr_base,
                         value_matrix[:, 0] + value_matrix[:, 2] + np.sqrt(np.power(value_matrix[:, 3], 2) + np.power(value_matrix[:, 5], 2)),
                         value_matrix[:, 0] + value_matrix[:, 2] - np.sqrt(np.power(value_matrix[:, 3], 2) + np.power(value_matrix[:, 5], 2)),
                         color='r',
                         alpha=.3)
    
        plt.legend(loc='best', ncol=1, title='')
        plt.axis([0.95, 2.05, 0, 4])
    
        plt.xticks(np.linspace(1, 2, 9), ['1.00', '', '1.25', '', '1.50', '', '1.75', '', '2.00'])
        plt.yticks(np.linspace(0, 4, 9))
    
        ax = plt.gca()
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='solid', alpha=0.5)
        ax.xaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    
        plt.xlabel("Ploidy")
        plt.ylabel("Abundance relative to total free molecule abundance in haploid")
    
        plt.savefig(os.path.join(Plotdir, "Abundance_vs_total_free_molecule_Size_%s.png" % key))
        plt.close()

def Plot_Osmotic_Pressure(y_min=8, y_max=18, yticks=8):

    ax = plt.gca()
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    ax.xaxis.grid(color='gray', linestyle='solid', alpha=0.5)

    plt.title('Cell diameter vs ploidy')
    plt.xlabel("Ploidy")
    plt.ylabel("Equivalent Diameter")

    plt.plot(arr_base,
             corrfactor * np.cbrt(means),
             '--k',
             label='',
             lw=2)

    plt.fill_between(arr_base,
                     corrfactor * np.cbrt(means - stds),
                     corrfactor * np.cbrt(means + stds),
                     color='k',
                     alpha=.3)


    plt.errorbar(ploidy_vs_size[sorting_index, 0],
                 ploidy_vs_size[sorting_index, 1]*corrfactor,
                 yerr=ploidy_vs_size[sorting_index, 2]*corrfactor,
                 fmt='ro--')

    plt.legend(loc='best', ncol=3)

    plt.savefig(os.path.join(Plotdir, "Cell_Diameter_vs_Ploidy.png"))
    plt.close()


    plt.suptitle(Suptitle)

    bplot1 = plt.boxplot([euploid_means, aneuploid_means],
                         patch_artist=True
                         )

    colors = ['lightblue', 'pink']
    for bplot in [bplot1]:
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

    vibration = 0.05*(2*np.random.rand(len(euploid_means))-1)
    plt.plot(1+vibration, euploid_means, 'ko')
    vibration = 0.05*(2*np.random.rand(len(aneuploid_means))-1)
    plt.plot(2+vibration, aneuploid_means, 'ko')

    plt.ylabel('Turgor pressure\nrelative to the euploid median')

    # add x-tick labels
    plt.xticks([y + 1 for y in range(2)], ['simulated\neuploid\nturgor',
                                           'simulated\naneuploid\nturgor'])

    #plt.show()
    plt.savefig(os.path.join(Plotdir, "Turgor_Pressure_euploid_aneuploid.png"))
    plt.close()
    #plot_abundances(pre_buckets)

def Plot_Sweep_Data(Sweepdatafile, sweep_parameter):
    """
    Plot the sweep data for the different sweep parameters
    """
    
    with open(Sweepdatafile) as fp: Sweep_Values = load(fp)
    arr_base = Sweep_Values["arr_base"]
    Sweep_Data = Sweep_Values["Sweep_Data"]
    Plot_Parameters = {
        "Title" : "Cell diameter vs ploidy vs %s" % sweep_parameter,
        "legend_title" : sweep_parameter,
        "base_label" : Base_Labels[sweep_parameter],
        "label_factor" : Label_Factor[sweep_parameter]
        }
    Plot_Parameters["Figure_File_Name"] = os.path.join(Plotdir, "%s.png" % Plot_Parameters["Title"].replace(" ", "_"))
    complex_size_swipe(Sweep_Data, arr_base, Plot_Parameters)


if __name__ == '__main__':

    Datadir = "data"
    Ploidyfile = os.path.join(Datadir, "ploidy_vs_size.dmp")
    ploidy_vs_size, corrfactor = Ploidy_Data(Ploidyfile)
    sorting_index = np.argsort(ploidy_vs_size[:, 0])
    Osmotic_Pressure_File = os.path.join(Datadir, "osmotic_pressure.dmp")

    Simdatafile = os.path.join(Datadir, "simulation_data.dmp")

    with open(Simdatafile) as fp: Sim_Data = load(fp)
    means = Sim_Data["means"]
    stds = Sim_Data["stds"]
    arr_base = Sim_Data["arr_base"]
    buckets = Sim_Data["buckets"]
    Observed_Volume, Predicted_Volume = Cell_Volumes(Sim_Data["arr_base"], Sim_Data["means"], Ploidyfile = Ploidyfile)
    press = Osmotic_Pressure(Observed_Volume, Predicted_Volume)
    aneuploid_means = press[1:-1]
    euploid_means = press[[0, -1], ]

    with open(Osmotic_Pressure_File) as fp:
        true_euploids_op, true_aneuploid_op = load(fp)

    normalizer_1 = np.median(euploid_means)
    normalizer_2 = np.median(true_euploids_op)

    true_euploids_op = np.array(true_euploids_op)/normalizer_2
    true_aneuploid_op = np.array(true_aneuploid_op)/normalizer_2

    euploid_means = euploid_means/normalizer_1
    aneuploid_means = aneuploid_means/normalizer_1
    Suptitle2 = "simulations medians ratio: 1:%.2f |" " true medians ratio: 1:%.2f" % (np.median(aneuploid_means)/np.median(euploid_means), np.median(true_aneuploid_op)/np.median(true_euploids_op))
    Suptitle = Suptitle1 + Suptitle2

    Plot_Osmotic_Pressure()
    plot_abundances(buckets)
    Sweepdatafile = os.path.join(Datadir, "Complex_size.dmp")
    Plot_Sweep_Data(Sweepdatafile, "complex size")

    #abundance_correlation_swipe()
    Sweepdatafile = os.path.join(Datadir, "Abundance_correlation.dmp")
    Plot_Sweep_Data(Sweepdatafile, "abundance correlation")

    #abundance_correlation_swipe()
    Sweepdatafile = os.path.join(Datadir, "Water_abundance.dmp")
    Plot_Sweep_Data(Sweepdatafile, "water abundance")

    Sweepdatafile = os.path.join(Datadir, "Ideality_correction.dmp")
    Plot_Sweep_Data(Sweepdatafile, "ideality correction")

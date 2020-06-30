#!/usr/bin/env python

import os
from random import shuffle
import numpy as np
from matplotlib import pyplot as plt
from pickle import load
import random
from matplotlib.cm import get_cmap
from statsmodels.nonparametric.smoothers_lowess import lowess
from calc_vol_press import Cell_Volumes, Osmotic_Pressure
from protein_abundance_preprocess import Ploidy_Data

Plotdir = "figures"

abundance_correlation = 0.7
large_complex_boost = 20
large_complex_correlation = 0.85

Suptitle1 = "large complex boost: x%.2f; large complex corr: %.2f, overall corr: %.2f\n" % (large_complex_boost, large_complex_correlation, abundance_correlation)


def diameter_plot_array(defautlt_params, base_array,
                        names_array=None, base_val=None, base_name=None,
                        array_vals_name='', y_min=8, y_max=15, yticks=8):
    """
    Generates a plot of states based on a variation of a single parameter.

    :param base_array: 1D numpy array
    :param names_array: 1D numpy string arrays

    :return:
    """
    plt.title('Cell diameter vs ploidy vs %s' % array_vals_name)

    base = np.linspace(0.0, 1.0, 20).tolist()
    arr_base = np.array(base)
    arr_base += 1

    ax = plt.gca()
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    ax.xaxis.grid(color='gray', linestyle='solid', alpha=0.5)

    plt.xlabel("Ploidy")
    plt.ylabel("Equivalent Diameter")

    cmap = get_cmap('brg')

    if names_array is None:
        names_array = base_array

    if base_val is not None and base_name is None:
        base_name = base_val

    # base array is assumed to be a 1D numpy array

    color_remap = np.linspace(0.5, 0.9, len(base_array))[np.argsort(base_array)[np.argsort(base_array)]]

    defautlt_params['base'] = base
    defautlt_params['total_partners'] = total_partners

    # for key, value in defautlt_params.iteritems():
    #     print key, value

    running_partial = core_sim_loop

    free_kwarg = ''
    for key, value in defautlt_params.iteritems():
        if value is not None:
            running_partial = partial(running_partial, **{key: value})
        else:
            free_kwarg = key

    # partial_function = partial(core_sim_loop, **defautlt_params)
    partial_function = running_partial

    # todo: form a partial function
    # => kwargs

    # reference plot
    if base_val is not None:

        means, stds, pre_buckets = partial_function(**{free_kwarg: base_val})

        plt.plot(arr_base,
                 corrfactor * np.cbrt(means),
                 '--k',
                 label=base_name,
                 lw=2)

        plt.fill_between(arr_base,
                         corrfactor * np.cbrt(means - stds),
                         corrfactor * np.cbrt(means + stds),
                         color='k',
                         alpha=.3)



    # the variation loop
    for i, (var_value, color, name) in enumerate(zip(base_array, color_remap, names_array)):

        # color generation from the 0.5-0.9 space
        color = cmap((color - 0.5) / (0.9 - 0.5) * 0.8 + 0.1)

        means, stds, pre_buckets = partial_function(**{free_kwarg: var_value})

        plt.plot(arr_base,
                 corrfactor * np.cbrt(means),
                 '--',
                 color=color,
                 label=name,
                 lw=2)

        plt.fill_between(arr_base,
                         corrfactor * np.cbrt(means - stds),
                         corrfactor * np.cbrt(means + stds),
                         color=color,
                         alpha=.3)

    # here what is plotted is varied as well.
    plt.legend(loc='best', ncol=3, title=array_vals_name)

    # and y-axis here needs to be adjusted dynamically as well
    #plt.axis([0.95, 2.05, y_min, y_max])
    plt.savefig(os.path.join(Plotdir, "Ploidy_vs_Cell_Diameter_%s.png" % array_vals_name.replace(" ", "_")))
    plt.close()

def complex_size_swipe():
    """
    Plots a distribution of complex sizes
    """
    Title = 'Cell diameter vs ploidy vs complex size'
    base = np.linspace(0.0, 1.0, 20).tolist()
    arr_base = np.array(base)
    arr_base += 1
    cmap = get_cmap('brg')
    plt.title(Title)
    plt.plot(arr_base, corrfactor*np.cbrt(arr_base), '--k', label='1', lw=2)
    
    ax = plt.gca()
    ax.set_axisbelow(True)
    
    ax.yaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    ax.xaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    
    plt.xticks(np.linspace(1, 2, 9), ['1.00', '', '1.25', '', '1.50', '', '1.75', '', '2.00'])
    plt.yticks(np.linspace(8, 15, 8), [8, '', 10, '', 12, '', 14, ''])
    
    plt.xlabel("Ploidy")
    plt.ylabel("Equivalent Diameter")
    
    complex_sizes = [2, 5, 10, 20, 40]
    
    for c_size in complex_sizes:
        c = cmap(c_size/40.*0.8+0.1)
    
        means, stds, buckets = core_sim_loop(base,[c_size], 0.75, 20, [])
    
        plt.plot(arr_base,
                 corrfactor*np.cbrt(means),
                 '--',
                 color=c,
                 label=c_size, lw=2)
    
        plt.fill_between(arr_base,
                         corrfactor*np.cbrt(means-stds),
                         corrfactor*np.cbrt(means+stds),
                         color=c,
                         alpha=.3)
    
    
    plt.legend(loc='best', ncol=3, title='complex size')
    #plt.axis([0.95, 2.05, 8, 15])
    plt.savefig(os.path.join(Plotdir, "%s.png" % Title.replace(" ", "_")))
    plt.close()

def plot_abundances(buckets):
    """
    Plots the abundance values of the proteins
    """
    for key, value_matrix in buckets.iteritems():
    
        plt.figure()
        plt.title('Molecules abundance vs ploidy % s for complex size %s' % (curr_corr, key))
    
        #print value_matrix
    
        norm = value_matrix[0, 0] + value_matrix[0, 2]
        value_matrix /= norm
    
        plt.plot(arr_base,
                 value_matrix[:, 0],
                 '-k',
                 label='complexes',
                 lw=2
                 )
    
        plt.plot(arr_base,
                 value_matrix[:, 2],
                 '-b',
                 label='free proteins',
                 lw=2
                 )
    
        plt.plot(arr_base,
                 value_matrix[:, 0] + value_matrix[:, 2],
                 '-r',
                 label='total molecules',
                 lw=2
                 )
    
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
    
        plt.savefig(os.path.join(Plotdir, "Abundance_vs_total_free_molecule.png"))
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

def abundance_correlation_swipe():
    """
    Runs the simulations for different complex abundances, alpha
    values and ideality correction factors
    """
    corr_vars = np.linspace(0.5, 0.9, 5)
    kwargs = {'abundance_correlation': None, 'repeats': 20, 'buckets': [5, 13], 'alpha': 1, 'ideality_correction': 1}
    diameter_plot_array(kwargs, corr_vars, corr_vars*100, 0, 0, 'complex abundance corr. (%)', 8, 15, 8)

    corr_vars = np.linspace(0.15, 0.95, 5)
    kwargs = {'abundance_correlation': 0.7, 'repeats': 20, 'buckets': [5, 13], 'alpha': None, 'ideality_correction': 1}
    diameter_plot_array(kwargs, corr_vars, corr_vars*100, 1.0, 100, 'alpha value. (%)', 4, 14, 11)
    
    corr_vars = np.linspace(0.5, 1.5, 6)
    kwargs = {'abundance_correlation': 0.7, 'repeats': 20, 'buckets': [5, 13], 'alpha': 1, 'ideality_correction': None}
    diameter_plot_array(kwargs, corr_vars, corr_vars*100, 1.0, 100, 'ideality correction factor. (%)', 8, 15, 8)


if __name__ == '__main__':

    Ploidyfile = "ploidy_vs_size.dmp"
    ploidy_vs_size, corrfactor = Ploidy_Data(Ploidyfile)
    sorting_index = np.argsort(ploidy_vs_size[:, 0])
    Osmotic_Pressure_File = "osmotic_pressure.dmp"

    Simdatafile = "simulation_data.dmp"

    with open(Simdatafile) as fp: Sim_Data = load(fp)
    means = Sim_Data["means"]
    stds = Sim_Data["stds"]
    arr_base = Sim_Data["arr_base"]
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
    complex_size_swipe()

    abundance_correlation_swipe()

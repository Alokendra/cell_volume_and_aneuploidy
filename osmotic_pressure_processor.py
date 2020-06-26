#!/usr/bin/env python

import os
from csv import reader as csv_reader
import numpy as np
from pickle import dump

type_index = []
master_table = []
Datadir = "data"
datafile = os.path.join(Datadir, "20171107_yeast.csv")


def clear_outlier(quadruplet, FDR=0.05, no_clear=True):

    potential_outliers = []

    if no_clear:
        return quadruplet

    for i in range(0, 3):
        potential_outlier = quadruplet[i]
        others_index = np.ones((4,), bool)
        others_index[i] = False
        all_others = quadruplet[others_index]
        lower, upper = t.interval(1-FDR, 2, np.mean(all_others), np.std(all_others))
        if potential_outlier > upper or potential_outlier < lower:
            potential_outliers.append(i)

    if len(potential_outliers) == 0:
        print 'no one cleared in %s' % quadruplet
        return quadruplet

    if len(potential_outliers) == 1:
        others_index = np.ones((4,), bool)
        others_index[potential_outliers[0]] = False
        print 'cleared %s in %s' % (quadruplet[potential_outliers[0]], quadruplet)
        return quadruplet[others_index]

    raise Exception('something went wrong with outlier calculation for the set %s; %s' %
                    (quadruplet, potential_outliers))


with open(datafile, 'rb') as source_file:
    reader = csv_reader(source_file)
    header = reader.next()
    for line in reader:
        strain_type = line[1][:1]
        strain_index = line[1][1:]
        type_index.append(strain_type)
        master_table.append((float(strain_index), float(line[2])))

type_index = np.array(type_index)
master_table = np.array(master_table)

aneuploid_selector = type_index == 'a'
euploid_selector = type_index == 'e'

euploid_raw = master_table[euploid_selector][:, 1]
aneuploid_raw = master_table[aneuploid_selector][:, 1]

aneuploid_clusters = np.unique(master_table[aneuploid_selector, 0])
euploid_clusters = np.unique(master_table[euploid_selector, 0])

aneuploid_cleared_clusters = []
euploid_cleared_clusters = []

for cluster in aneuploid_clusters:
    aneuploid_cleared_clusters.append(clear_outlier(master_table[np.logical_and(master_table[:,
                                                                         0] == cluster,
                                               aneuploid_selector), 1]))

for cluster in euploid_clusters:
    euploid_cleared_clusters.append(clear_outlier(master_table[np.logical_and(master_table[:,
                                                                        0] == cluster,
                                             euploid_selector), 1]))

euploid_means = np.mean(euploid_cleared_clusters, axis = 1).tolist()
aneuploid_means = np.mean(aneuploid_cleared_clusters, axis = 1).tolist()

print "means ratio: 1:%f" % (np.mean(aneuploid_means)/np.mean(euploid_means))

Cluster_Data = {"euploid_cleared_clusters" : euploid_cleared_clusters, "aneuploid_cleared_clusters" : aneuploid_cleared_clusters}
Cluster_Mean_Data = { "euploid_means" : euploid_means, "aneuploid_means" : aneuploid_means }
Cluster_Raw_Data = { "euploid_raw" : euploid_raw, "aneuploid_raw" : aneuploid_raw }

with open('osmotic_pressure.dmp', 'w') as fp:
    dump((euploid_means, aneuploid_means), fp)
with open('osmotic_pressure_all.dmp', 'w') as fp:
    dump(Cluster_Data, fp)

from osmotic_pressure_plots import smooth_histogram, plot_ploidy_means, boxplots

smooth_histogram(Cluster_Raw_Data, Figfilename = "Osmotic_Pressure_Histogam.png")
plot_ploidy_means(Cluster_Data, Figfilename = "Osmotic_Pressure_Mean.png")
boxplots(Cluster_Mean_Data, Figfilename = "Osmotic_Pressure_Box_Plots")
# print euploids_raw

# smooth_histogram(euploids_raw)
# smooth_histogram(aneuploids_raw, 'r')
#
# plt.show()



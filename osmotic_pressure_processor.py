#!/usr/bin/env python

import os
from csv import reader as csv_reader
import numpy as np
from pickle import dump


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

def Build_Table(datafile):
    """
    Build master table
    """
    type_index = []
    master_table = []
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

    return type_index, master_table


def Collect_Data(type_index, master_table):
    """
    Collects the data for euploids and aneuploids
    """
    Ploidy_Data = {"Euploids" : {}, "Aneuploids" : {}}
    Selector = {"Euploids" : "a", "Aneuploids" : "e"}

    for key in Ploidy_Data:
        Ploidy_Data[key]["selector"] = type_index == Selector[key]
        Ploidy_Data[key]["raw"] = master_table[Ploidy_Data[key]["selector"]][:, 1]
        Ploidy_Data[key]["clusters"] = np.unique(master_table[Ploidy_Data[key]["selector"], 0])
        Cleared_Clusters = []
        for cluster in Ploidy_Data[key]["clusters"]:
            Cleared_Clusters.append(clear_outlier(master_table[np.logical_and(master_table[:, 0] == cluster, Ploidy_Data[key]["selector"]), 1]))
        Ploidy_Data[key]["cleared_clusters"] = Cleared_Clusters
        Ploidy_Data[key]["cluster_means"] = np.mean(Cleared_Clusters, axis = 1).tolist()

    return Ploidy_Data

if __name__ == '__main__':

    Datadir = "data"
    datafile = os.path.join(Datadir, "20171107_yeast.csv")

    Osmotic_Pressure_Dumpfile = os.path.join(Datadir, "osmotic_pressure.dmp")
    Osmotic_Pressure_AllDumpfile = os.path.join(Datadir, "osmotic_pressure_all.dmp")
    type_index, master_table = Build_Table(datafile)
    Ploidy_Data = Collect_Data(type_index, master_table)

    print "means ratio: 1:%f" % (np.mean(Ploidy_Data["Euploids"]["cluster_means"])/np.mean(Ploidy_Data["Aneuploids"]["cluster_means"]))


    with open(Osmotic_Pressure_Dumpfile, 'w') as fp:
        dump((Ploidy_Data["Euploids"]["cluster_means"], Ploidy_Data["Aneuploids"]["cluster_means"]), fp)
    with open(Osmotic_Pressure_AllDumpfile, 'w') as fp:
        dump(Ploidy_Data, fp)


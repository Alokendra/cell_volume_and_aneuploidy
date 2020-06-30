#!/usr/bin/env python

# This script creates the plots from the osmotic pressure data

import os, pickle
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde, t, sem

Textlabel = lambda strng : strng.replace("_", " ").capitalize()

def Saveorshow(plt, Figfilename):
    """
    Save or show figure
    """
    if Figfilename:
        plt.savefig(Figfilename)
        plt.close()
    else:
        plt.show()

def smooth_histogram(Cluster_Raw, Figfilename = ""):
    """
    Plots the histograms from raw data
    """
    colorlist = ['k','r','b']
    Labels = ["euploid_raw","aneuploid_raw"]
    for ii, key in enumerate(Labels):
        data = Cluster_Raw[key]
        fltr = np.logical_not(np.isnan(data))
        density = gaussian_kde(data[fltr].flatten())
        xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
        plt.plot(xs, density(xs), colorlist[ii])

    plt.legend(map(Textlabel, Labels))

    Saveorshow(plt, os.path.join(Plotdir, Figfilename))



def plot_ploidy_means(Cluster_Data, Figfilename = ""):
    """
    Creates the mean clusters for aneuploid and euploids
    from external data file
    """
    Plots = []
    Fmt = dict(zip(Cluster_Data.keys(), ["rs--","ks--"]))

    for key, value in Cluster_Data.items():
        for i, cluster in enumerate(value):
            line, upper, lower = plt.errorbar(i, np.mean(cluster), yerr=sem(cluster), fmt=Fmt[key])
        line.set_label(Textlabel(key))

    plt.legend()

    Saveorshow(plt, os.path.join(Plotdir, Figfilename))


def boxplots(Cluster_Mean_Data, Figfilename = ""):
    """
    Gets the box plot of data
    """
    Labels = ["euploid_means", "aneuploid_means"]
    Markers = ["ko", "ro", "bo"]
    for ii, key in enumerate(Labels):
        vibration = 0.05*(2*np.random.rand(len(Cluster_Mean_Data[key]))-1)
        plt.plot(ii+1+vibration, Cluster_Mean_Data[key], Markers[ii])

    plt.legend(map(Textlabel, Labels))
    plt.boxplot([Cluster_Mean_Data["euploid_means"], Cluster_Mean_Data["aneuploid_means"]])

    Saveorshow(plt, os.path.join(Plotdir, Figfilename))

if __name__ == '__main__':

    Datadir = "data"
    Dumpfile = os.path.join(Datadir, "osmotic_pressure_all.dmp")
    Plotdir = "figures"

    with open(Dumpfile) as fp: Ploidy_Data = pickle.load(fp)
    
    Cluster_Data = {"euploid_cleared_clusters" : Ploidy_Data["Euploids"]["cleared_clusters"], "aneuploid_cleared_clusters" : Ploidy_Data["Aneuploids"]["cleared_clusters"]}
    Cluster_Mean_Data = { "euploid_means" : Ploidy_Data["Euploids"]["cluster_means"], "aneuploid_means" : Ploidy_Data["Aneuploids"]["cluster_means"]}
    Cluster_Raw_Data = { "euploid_raw" : Ploidy_Data["Euploids"]["raw"], "aneuploid_raw" : Ploidy_Data["Aneuploids"]["raw"]}


    smooth_histogram(Cluster_Raw_Data, Figfilename = "Osmotic_Pressure_Histogam.png")
    plot_ploidy_means(Cluster_Data, Figfilename = "Osmotic_Pressure_Mean.png")
    boxplots(Cluster_Mean_Data, Figfilename = "Osmotic_Pressure_Box_Plots")

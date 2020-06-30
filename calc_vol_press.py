#!/usr/bin/env python

import os
from random import shuffle
import numpy as np
from pickle import load, dump
import random
from scipy.interpolate import interp1d
from protein_abundance_preprocess import Ploidy_Data

Parameters = {
    "Vm" : 18.03e-3,
    "Pi" : 0.05e6,
    "T" : 293
    }

Abundance_Factor = lambda Pi, Vm, T = 293: np.exp(-Pi*Vm/(8.314*T))

def Cell_Volumes(arr_base, means, Ploidyfile = "ploidy_vs_size.dmp"):
    """
    Calculates the osmotic pressure using the protein
    abundances
    """
    ploidy_vs_size, corrfactor = Ploidy_Data(Ploidyfile)
    interp_prediction = interp1d(arr_base, corrfactor*np.cbrt(means), kind='cubic')

    sorting_index = np.argsort(ploidy_vs_size[:, 0])

    Observed_Volume = (ploidy_vs_size[sorting_index, 1]*corrfactor)**3
    Predicted_Volume = (interp_prediction(ploidy_vs_size[sorting_index, 0]))**3

    return Observed_Volume, Predicted_Volume

def Osmotic_Pressure(Observed_Volume, Predicted_Volume):
    """
    Calculates the osmotic pressure from predicted and
    observed volumes
    """
    alpha_0 = Abundance_Factor(Parameters["Pi"], Parameters["Vm"], Parameters["T"])
    undervolume = Observed_Volume/Predicted_Volume

    press = []
    for value in undervolume:
        alpha = 1 - (1 - alpha_0)/value
        pressure = -8.314*Parameters["T"]/Parameters["Vm"]*np.log(alpha)
        # print value, alpha, pressure, pressure/0.05e6
        if not np.isnan(pressure):
            press.append(pressure)

    press = np.array(press)

    return press

if __name__ == "__main__":

    Ploidyfile = "ploidy_vs_size.dmp"
    ploidy_vs_size, corrfactor = Ploidy_Data(Ploidyfile)

    Simdatafile = "simulation_data.dmp"

    with open(Simdatafile) as fp: Sim_Data = load(fp)
    Observed_Volume, Predicted_Volume = Cell_Volume(Sim_Data["arr_base"], Sim_Data["means"], corrfactor)
    press = Osmotic_Pressure(Observed_Volume, Predicted_Volume)

    with open("volume_pressure_data.dmp", "w") as fp:
        dump({ "Observed_Volume" : Observed_Volume,
               "Predicted_Volume" : Predicted_Volume,
               "press" : press }, fp)

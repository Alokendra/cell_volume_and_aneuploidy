---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.5.0
  kernelspec:
    display_name: Python (py2dev)
    language: python
    name: py2dev
---

<!-- #region slideshow={"slide_type": "slide"} -->
## Cell Volume dependence on protein abundance and osmotic stress

This document is created to explore some of the results from Tsai et. al. biophysical model which that derives a relationship between cell volume and ploidy by showing how changes in the *protein abundances* and *composition* in aneuploid cells creates a hypo osmotic stress which causes a change in cell volume. Here we show the main equations of the model and its application to show how cell volume changes as a function of ploidy under different combinations of biophysical parameters of the cell. For the ease of exploration sliders are created for some of the key parameters of the model that can be changed interactively.
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
### Model Equations

Considering transport of water and by balancing the chemical potential of water inside and outside the cell one gets

$$\mu_{H_2O}^{in} - RT \ln c_{H_2O}^{in} - \int_P^{P + \Pi}V_mdp = \mu_{H_2O}^{out} + \Delta\mu_{H_2O}^{memb} - RT \ln c_{H_2O}^{out}$$
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
At equilibrium, by equating the chemical potential of water inside and outside one can write

$$\Pi = -\frac{RT}{V_m}\ln{\left(\frac{c_{H_2O}^{in}}{c_{H_2O}^{out}c_{H_2O}^{memb}}\right)} = -\frac{RT}{V_m}\ln{\alpha}$$

In the above equation, $\alpha$ is the water abundance factor.
Considering that water is the dominant component of the cell the effective volume of the cell is then proportional to number of proteins in the cell times the water abundance factor $\alpha$

$$V_{eff} = V - V_m \propto N_p\alpha$$
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
If there are $N_k$ proteins of type k (where $k = 1,2,\ldots,N$) in the cell and $N_c$ complexes of type c (where $c = 1,2,\ldots,M$ then the total number of free proteins is given by [Tsai 2019][tsai2019hypo]

$$N_p = \sum_{k=1}^{N}N_k - \sum_{c=1}^{M}\min_{k \in S_c}(N_k)N_c$$
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
Here, the second term considers that the number of complexes of type $c$ is equal to the number of its least abundant component.
This holds if every protein participates in a single complex and every complex assembles fully which are assumed here.
Now if we have aneuploid cells then for particular proteins their amount may double due to duplication of gene coding.
If $A_k$ is the abundance of protein $k$ and $P_k$ is the probability of gene duplication then the number of proteins can be written as $N_k = A_kP_k$. Then the above equatio becomes

$$V_{eff} = V - V_m = \alpha \left(\sum_{k=1}^{N}A_kP_k - \sum_{c=1}^M\min_{k \in S_c}(A_kP_k)N_c\right)$$
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
To obtain the effective volume one needs estimates of relative abundances of proteins and average sizes of protein complexes. In [Tsai 2019][tsai2019hypo] the relative abundances were obtained from PAX-DB project dataset and average protein complex size was estimated using two methods

- Average number of protein-protein interactions.
- Measurements of number of complex size distributions.

A correction factor was also introduced for protein abundance to account for the correlation of abundances in proteins that form complex.


[tsai2019hypo]: #References
<!-- #endregion -->

<!-- #region slideshow={"slide_type": "slide"} -->
### Simulation

The simulations were conducted for different sets of parameters. First the important modules were imported. Numpy and matplotlib needs to be installed.
<!-- #endregion -->

```python slideshow={"slide_type": "skip"}
import os, sys
from random import shuffle
import numpy as np
from pickle import load, dump
from combined_simulations import Get_Abundance_Data, core_sim_loop, average_simulation_data
from protein_abundance_preprocess import Ploidy_Data
from protein_abundances import generate_complex_ids, last_complex_adjustment, align_complex_abundances,\
sorted_complex_abundances
from simulation_plots import Plot_Sweep_Data
```

<!-- #region slideshow={"slide_type": "slide"} -->
#### Complex size distribution

The first case we consider is the simplest, where all complexes are present at the same abundance (6). The different parameters of the models are set below
<!-- #endregion -->

```python slideshow={"slide_type": "subslide"}
Sweep_Parameters = {
    "complex size" : {
        "sweep_elem" : "total_partners",
        "total_partners" : [6],
        "abundance_correlation" : 0.7,
        "alpha" : 1,
        "ideality_correction" : 1
        }
}
```

<!-- #region slideshow={"slide_type": "skip"} -->
We can initialize some other variables as below
<!-- #endregion -->

```python slideshow={"slide_type": "skip"}
Datadir = "data"
Plotdir = "figures"
Datastatfile = os.path.join(Datadir, "data_stats_dump.dmp")
Interaction = "Paxdb"
abundance_range, total_partners = Get_Abundance_Data(Datastatfile, Interaction = Interaction)

base = np.linspace(0.0, 1.0, 20).tolist()
arr_base = np.array(base) + 1
```

<!-- #region slideshow={"slide_type": "skip"} -->
Then we can create the different cases as below
<!-- #endregion -->

```python slideshow={"slide_type": "skip"}
def sweep_parameter(base, sweep_name):
    """
    Performs a sweep over complex sizes
    """
    sweep_key = Sweep_Parameters[sweep_name]["sweep_elem"]
    sweep_values = Sweep_Parameters[sweep_name][sweep_key]
    Data = { key : [] for key in sweep_values }
    arr_base = np.array(base) + 1
    Sweepfile = os.path.join("data", "%s.dmp" % sweep_name.capitalize().replace(" ","_"))

    for key in sweep_values:

        Sweep_Parameters[sweep_name][sweep_key] = key
        complex_contents = generate_complex_ids([Sweep_Parameters[sweep_name]["total_partners"]], len(abundance_range))
        complex_contents = last_complex_adjustment(complex_contents, len(abundance_range))
        aligned_abundances = align_complex_abundances(complex_contents, abundance_range, abundance_correlation = Sweep_Parameters[sweep_name]["abundance_correlation"])
        aligned_abundances = sorted_complex_abundances(aligned_abundances, complex_contents[-1], abundance_range, abundance_correlation = Sweep_Parameters[sweep_name]["abundance_correlation"])
        re_runs, buckets = core_sim_loop(base, complex_contents, aligned_abundances)
        means, stds, pre_buckets = average_simulation_data(re_runs, buckets, alpha = Sweep_Parameters[sweep_name]["alpha"], ideality_correction = Sweep_Parameters[sweep_name]["ideality_correction"])
        Data[key] = (means, stds, pre_buckets)

    with open(Sweepfile, "w") as fp: dump({"Sweep_Data" : Data, "arr_base" : arr_base}, fp)
```

<!-- #region slideshow={"slide_type": "slide"} -->
Next we can run the sweep over the complex sizes (in this case there is only one, 6)
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
sweep_parameter(base, "complex size")
```

<!-- #region slideshow={"slide_type": "slide"} -->
Then we can import the plotting module and get the plots as below
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
Sweepdatafile = os.path.join(Datadir, "Complex_size.dmp")
from IPython.display import Image, display
from ipywidgets import interactive
import ipywidgets as widgets

def plot_complex_size(Size=2):
    Sweep_Parameters["complex size"]["total_partners"] = [Size]
    sweep_parameter(base, "complex size")
    Plot_Sweep_Data(Sweepdatafile, "complex size")
    display(Image(filename='figures/Cell_diameter_vs_ploidy_vs_complex_size.png'))

p = interactive(plot_complex_size, Size = widgets.IntSlider(min=5, max=50, step=5, value=5))
p
```

<!-- #region slideshow={"slide_type": "slide"} -->
Next we will change the complex numbers over a range instead of a single value and rerun the simulation again
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
def plot_complex_sweep(abundance_correlation, alpha):
    Sweep_Parameters["complex size"]["total_partners"] = [2, 5, 10, 20, 40]
    Sweep_Parameters["complex size"]["abundance_correlation"] = abundance_correlation
    Sweep_Parameters["complex size"]["alpha"] = alpha
    sweep_parameter(base, "complex size")
    Sweepdatafile = os.path.join(Datadir, "Complex_size.dmp")
    Plot_Sweep_Data(Sweepdatafile, "complex size")
    display(Image(filename='figures/Cell_diameter_vs_ploidy_vs_complex_size.png'))
    
q = interactive(plot_complex_sweep, 
                abundance_correlation = widgets.FloatSlider(min=0.5, max=0.9, step=0.1, value=0.7),
                alpha = widgets.FloatSlider(min=0.2, max=1.0, step=0.2, value=1.0))
q
```

<!-- #region slideshow={"slide_type": "slide"} -->
#### Abundance Correlation Factor

In a similar way we can set the the abundance correlation factor as below
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
Sweep_Parameters["abundance correlation"] = {
        "sweep_elem" : "abundance_correlation",
        "total_partners" : 2,
        "abundance_correlation" : np.linspace(0.5, 0.9, 5).tolist(),
        "alpha" : 1,
        "ideality_correction" : 1
        }
```

<!-- #region slideshow={"slide_type": "slide"} -->
Now we can run the simulation and plot the results in a similar way
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
def plot_abundance_correlation(total_partners, alpha):
    Sweep_Parameters["abundance correlation"]["total_partners"] = total_partners
    Sweep_Parameters["abundance correlation"]["abundance_correlation"] = np.linspace(0.5, 0.9, 5).tolist()
    Sweep_Parameters["abundance correlation"]["alpha"] = alpha
    sweep_parameter(base, "abundance correlation")
    Sweepdatafile = os.path.join(Datadir, "Abundance_correlation.dmp")
    Plot_Sweep_Data(Sweepdatafile, "abundance correlation")
    display(Image(filename='figures/Cell_diameter_vs_ploidy_vs_abundance_correlation.png'))
    
r = interactive(plot_abundance_correlation, 
                total_partners = widgets.IntSlider(min=5, max=30, step=5, value=5),
                alpha = widgets.FloatSlider(min=0.2, max=1.0, step=0.2, value=1.0))
r
```

<!-- #region slideshow={"slide_type": "slide"} -->
#### Water abundance factor

Next we can change the water abundance factor and see its effect
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
Sweep_Parameters["water abundance"] = {
        "sweep_elem" : "alpha",
        "total_partners" : 5,
        "abundance_correlation" : 0.5,
        "alpha" : np.linspace(0.15, 1.0, 6).tolist(),
        "ideality_correction" : 1
        }
```

```python slideshow={"slide_type": "slide"}
def plot_water_abundance(total_partners, abundance_correlation):
    Sweep_Parameters["water abundance"]["total_partners"] = total_partners
    Sweep_Parameters["water abundance"]["abundance_correlation"] = abundance_correlation
    Sweep_Parameters["water abundance"]["alpha"] = np.linspace(0.15, 1.0, 6).tolist()
    sweep_parameter(base, "water abundance")
    Sweepdatafile = os.path.join(Datadir, "Water_abundance.dmp")
    Plot_Sweep_Data(Sweepdatafile, "water abundance")
    display(Image(filename='figures/Cell_diameter_vs_ploidy_vs_water_abundance.png'))
    
s = interactive(plot_water_abundance, 
                total_partners = widgets.IntSlider(min=5, max=30, step=5, value=5),
                abundance_correlation = widgets.FloatSlider(min=0.5, max=0.9, step=0.1, value=0.7))
s
```

<!-- #region slideshow={"slide_type": "slide"} -->
#### Ideality Correction Factor

Finally we can change the ideality correction factor as below
<!-- #endregion -->

```python slideshow={"slide_type": "slide"}
Sweep_Parameters["ideality correction"] = {
        "sweep_elem" : "ideality_correction",
        "total_partners" : 20,
        "abundance_correlation" : 0.7,
        "alpha" : 1,
        "ideality_correction" : np.linspace(0.5, 1.5, 6).tolist()
        }
```

```python slideshow={"slide_type": "slide"}
def plot_ideality_correction(total_partners, abundance_correlation, alpha):
    Sweep_Parameters["ideality correction"]["total_partners"] = total_partners
    Sweep_Parameters["ideality correction"]["abundance_correlation"] = abundance_correlation
    Sweep_Parameters["ideality correction"]["alpha"] = alpha
    Sweep_Parameters["ideality correction"]["ideality_correction"] = np.linspace(0.5, 1.5, 6).tolist()
    sweep_parameter(base, "ideality correction")
    Sweepdatafile = os.path.join(Datadir, "Ideality_correction.dmp")
    Plot_Sweep_Data(Sweepdatafile, "ideality correction")
    display(Image(filename='figures/Cell_diameter_vs_ploidy_vs_ideality_correction.png'))
    
t = interactive(plot_ideality_correction, 
                total_partners = widgets.IntSlider(min=5, max=30, step=5, value=5),
                abundance_correlation = widgets.FloatSlider(min=0.5, max=0.9, step=0.1, value=0.7),
                alpha = widgets.FloatSlider(min=0.2, max=1.0, step=0.2, value=1.0))
t
```

<!-- #region slideshow={"slide_type": "slide"} -->
### References
 
1. Hung-Ji Tsai, Anjali R Nelliat, Mohammad Ikbal Choudhury, Andrei Kuchar-
avy, William D Bradford, Malcolm E Cook, Jisoo Kim, Devin B Mair, Sean X
Sun, Michael C Schatz, and Rong Li. Hypo-osmotic-like stress underlies general
cellular defects of aneuploidy. Nature, 570(7759):117–121, 2019.
<!-- #endregion -->

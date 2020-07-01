# Cell Volume and Aneuploidy

This repository is a modified version of the model originally
published in "Tsai, Hung-Ji, et al. Nature 570.7759 (2019):
117-121.". This is a biophysical model of cell volume and osmotic
pressure dependence on the ploidy of a cell. This repository
provides an interactive front-end to the model through a [Jupyter
Notebook](https://jupyter.org).

## Cloning

This repository data files are maintained using [git-lfs](https://github.com/git-lfs/git-lfs). In order to clone this from command line `git-lfs` needs to be installed in the target system. The latest binary for linux can be obtained from [here](https://github.com/git-lfs/git-lfs/releases/download/v2.11.0/git-lfs-linux-amd64-v2.11.0.tar.gz). After unpacking, simply run the `install.sh` script provided. After installation the repository can be cloned normally.


## Installation

To install and view the jupyter notebook create a copy of the
this repository and run the helper shell script `Virtualenv.sh`
as below

```sh
$ chmod +x Virtualenv.sh
$ ./Virtualenv.sh
```

This script will setup a virtual environment, install the
necessary packages in it and (optionally) open the notebook file
in a browser window from where it can be run.

NOTE: By default the notebook is set as untrusted. It is
recommended to set it as trusted before running the code.

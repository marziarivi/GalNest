# GalNest
[Radio Galaxy Detection in the Visibility Domain](https://arxiv.org/abs/1810.12930), A. Malyali, M. Rivi, F.B. Abdalla, J.D. McEwen, 2019, MNRAS 486(2) 2695–2704

Bayesian model-fitting method working in the Fourier domain adopting a single SF galaxy model (the exponential profile, i.e. Sersic of index 1). The resulting multimodal posterior distribution is sampled using a [multimodal nested sampling](https://ccpforge.cse.rl.ac.uk/gf/project/multinest/) algorithm. 

This is the C version of the [python code](https://github.com/amalyali/RadioGalFit/tree/master/GalNest). MPI Parallelization MPI can be enabled in the Makefile (each MPI task will read a different spectral window of the dataset that must be split in independent MS files).

MultiNest output files must be processed to remove fake modes using the python code _Modes_selection.py_, specifying the flux/arcsec threshold.

## Installation

### GalNest dependencies
1. [casacore](https://github.com/casacore/casacore) library for reading Measurement Set file format.
2. [MultiNest](https://github.com/JohannesBuchner/MultiNest) library.  
3. MPI library for multi-node parallelization

Update the Makefile with the path to the dependencies.

### Modes_selection dependencies
[pyMultiNest](https://github.com/JohannesBuchner/PyMultiNest) package for reading MultiNest output files.

# GimMIRI

An Nvidia/CUDA C based simulator that models the internal diffraction within detector substrates, which result in a cross-like feature around bright sources and a broad halo. This effect has been observed for Si:As type blocked impurity band infrared as well as regular CCD detectors. Our code is a Quantum-Electrodynamics based model and does not therefore assume any reflection laws or physics, other than the Fresnel laws of reflection at the backside of the detectors. The model therefore provides a smooth transition between far and near-field approximations. The downside of the QED model is that it is numerically expensive, possibly taking many days for a single calculation to execute. The code, therefore, is GPU based.

# Installing

The github package contains the source code, test files, example PBS files for executing on a computer cluster, and a cfitsio directory. 

Once downloaded, cfitsio needs to be compiled first, unless the system already has cfitsio installed.

$ cd cfitsio

$ ./configure

$ make

$ cd ../

If a system-wide cfitsio is used, the Makefile in src needs to be edited to point to proper libraries and include locations. Otherwise, a simple make command will compile all the necessary components of GimMIRI.

$ make

# Running GimMIRI

GimMIRI takes an input parameter file that contains all the information it needs to run. Example parameter files can be found in the params directory. The variables are explained in the parameter file and all file locations are given relative to the DiskDyn directory.

$ ./GimMIRI -i params/Reference.param

A single example MIRI PSF is placed in the PSFlib directory, with a scaling factor of one. It is recommended to use much larger (oversized) PSF for higher accuracy, which can be generated using WebbPSF. 

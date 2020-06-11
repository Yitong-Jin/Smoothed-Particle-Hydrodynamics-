# ACSE-4-SPH

[Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) (SPH) is a meshless
method for solving the Navier-Stokes equation, in which fluid properties are stored on Lagrangian fluid particles (i.e. on
particles which move with the fluid flow). The particles interact to generate values across the entire fluid domain through
continuous smoothing kernels. 

As the SPH method is meshless and Lagrangian, it is ideal for solving problems involving fluid flow with interfaces and free 
surfaces. This tool implements the SPH method in C++ to solve wave generation in a lock-release/dam-break problem.

### Compilation/Installation Guide

This model uses C++ with OpenMP. 

To import the code (when in the chosen directory) go into the terminal and enter:
```
git clone git@github.com:acse-2019/acse-4-sph-humphrey.git
```
To check whether the ouput files are correct, Python must also be installed with the vtk library. To install, write the following into the command line:
```
pip install vtk
```
To read the output files, use Paraview. Download Paraview (v 4.0.1) at https://www.paraview.org/download/.

### User instructions

To use forward Euler

To use improved Euler

To use flat domain

To use beach domain

### Documentation

The code includes deoxygen documentation. On systems with deoxygen installed, this can be built by create a configuration file and running deoxygen ?????

### Testing

The tool includes tests, which you can use to check its operation on your system. With the code compiled, these can be run 
with

```
python -m pytest SPH_Snippet.cpp ???
```

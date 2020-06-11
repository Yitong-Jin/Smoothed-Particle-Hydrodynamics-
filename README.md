# ACSE-4-SPH Team Humphrey

![](wave_velocity.gif)

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

By default, our simulations use the Predictor-Corrector timestepping method. If you wish to use Forward Euler, there is a nested repository to use this timestepping method. The scripts are executed in the same manner from the base directory, except including the additional `Forward_Euler_repo` in the directory path.

The model allows the user to run on box containers and containers that replicate beach environments to model shoaling waves. This is selected by the user at the start of the `SPH_Snippet.cpp` script, where they can select `beach` or `no beach` in the global environment by uncommenting appropriately. The user can then define how long the model runs for by changing the `max_time` variable, and change how often a file is written using the `print_every` variable (as the timestep is adaptive, the exact printout interval is not exact but is as close to constant intervals as possible.

The user can then visualize the resulting VTP files in ParaView, or use the Python post-processing script `generate_animation.py` to generate a **.mp4** file to be visualized on your media viewer of choice.

### Running file

* Windows (MVSC) Users:
  - Build and run within main folder, changing the directory structure as prescribed below (Linux/iOS Users).
  
 * Linux/iOS Users:
  - Compile with g++ compiler from GNU from base directory:
  Improved Euler:
  `g++ -std=c++14 -fopenmp -o SPH ./includes/*.cpp`
  Foward Euler
  `g++ -std=c++14 -fopenmp -o SPH ./Forward_Euler_repo/includes/*.cpp` 
  Run using the specified target (from `-o SPH`):
  
  `./SPH`
  
  The output files are then accessed in `./tests/` or `./Forward_Euler_repo/tests/` for Predictor-Corrector (Improved Euler) and Forward Euler, respectively.

### Documentation

On systems with deoxygen installed, this can be built by running doxygen and choosing a file that contains doxygen documentation. A folder named "Documentation" holds the html and latex files that was formed by with doxygen and SPH_2D.h. This provides an easy way to read what each class and function does.

The folder also contains SPH_2D.h with doxygen documentation, as well as SPH_2D.cpp with doxygen documentation. Using doxygen on the .cpp is not recommended, however, unless the user wants to read through the steps (in words) for each function. Reading the comments along with the code without doxygen would be easier.

Further information can be found in Doxyfile - it describes the settings to be used by the documentation system.


### Testing

The tool includes test, which you can use to check correctness of output files with running "python_tests.py". To use the test, enter "tests/" directory and it can be run with:

```
python python_tests.py
```

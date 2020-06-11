#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

class SPH_main;

class SPH_particle
{
public:

	SPH_particle();

	double x[2][2], v[2][2];			//position and velocity
	double rho[2], P;					//density and pressure
	double mass;						// mass in kg
	double rho0;						// density in kg / m^3
	double dv[2];						// acceleration in x and y directions in m / s^2
	double drho;						// (big) D - rate of density change
	double visc;						// viscosity in Pa s
	double gamma;						// Stiffness
	bool boundary;						// indicates whether particle is part of boundary or not
	int part_num;						// number of particle in sub_grid
	double rho_new;						// temporary rho storage used during smoothing
	int step;							// step incrementor of predictor-corrector method

	static SPH_main* main_data;			//link to SPH_main class so that it can be used in calc_index

	int list_num[2];					//index in neighbour finding array

	// Function to calculate the row and col index of the particle's sub-grid
	void calc_index(void);

	// Function to calculate B, where B is the reference density of water
	double get_B(void);

	// Function to get P, where P is pressure - calculated using the Tait Equation
	double get_P(double& rho);

	// Function to get the differential (dwdr) of the smoothing kernel
	double get_dwdr(double& r, double& h);

	// Function that calculates the acceleration of a particle based on summation using the smoothing kernel's gradient
	// Will be used in the neighbour_iterate_stencil
	// Also updates dt_F
	void NS_acceleration_simultaneous(SPH_particle* other_particle, double* e_ij, double* v_ij, double& dist);

	// Function that calculates the acceleration of a particle based on summation using the smoothing kernel's gradient
	// Will be used in the neighbour_iterate
	// Also updates dt_F
	void NS_acceleration(SPH_particle* other_particle, double* e_ij, double* v_ij, double dist);

	// Function that calculates the density change of a particle over time based on summation using the smoothing kernel's gradient
	// Will be used in the neighbour_iterate_stencil
	// Also updates dt_A
	void NS_density_simultaneous(SPH_particle* other_particle, double* e_ij, double* v_ij, double& dist);

	// Function that calculates the density change of a particle over time based on summation using the smoothing kernel's gradient
	// Will be used in the neighbour_iterate
	// Also updates dt_A
	void NS_density(SPH_particle* other_particle, double* e_ij, double* v_ij, double dist);

	// W is the smoothing kernel - has compact support and it is chosen to be normal in this SPH model
	double get_W(double& r, double& h);

	// Function to calculate the repulsive force 
	double repulsion_func(double& d, double& n);	// Defines a function that's repulsive at the boundaries and forceless elsewhere

	//Function to apply the repulsive forces onto the fluid particles
	void apply_boundary_forces(string keyword);

	// Function to update the properties of the particles after every iteration using Forward Euler method
	void update(string keyword);					// Forward

	// Function to update the properties of the particles after every iteration using Improved Euler method
	void update_pc(int step, string keyword);		// Predictor-Corrector

	// Set half-step increments (required by update_pc)
	void set_step(int step_);
};


class SPH_main
{
public:
	SPH_main();

	// Function that set values for the entire space
	void set_values(void);

	// Function to initialise a grid onto the space
	void initialise_grid(void);

	// Function to place both the fluid and boundary particles into the grid
	void place_points(double* min, double* max);

	// Function to place both the fluid and boundary particles into the grid for a beach domain
	void place_points_beach(double* min, double* max);

	// Function that allocates all the points to the search grid (assumes that index has been appropriately updated)
	void allocate_to_grid(void);

	// Function to calculate the acceleration and density by iterating over particles within 2h of the chosen particle
	void neighbour_iterate(SPH_particle* part);

	// Function that smooths densities of all the particles
	void smoothing();

	// Function to smooth the density for a particular particle
	// Also updates dt_CFL
	void smoothing_neighbour_itr(SPH_particle* part);

	// Function to calculate the acceleration and density by iterating over particles within 2h of the chosen particle using a stencil
	// Also updates dt_CFL
	void neighbour_iterate_stencil(SPH_particle* part);

	// Function to update the time-step according to the properties of the current time
	void set_adaptive_dt();

	double h_fac;							// set in set values - links h to dx^2 (1.3)

	double dx;								//particle initial spacing

	double g;								// acceleration due to gravity

	double c_0;								// speed of sound

	double dt;								// set in set_values
	double dt_CFL;							// Courant-Freidrichs-Lewy Criterion
	double dt_F;							// Force criterion
	double dt_A;							// Acoustic criterion

	double h;								//smoothing length

	double min_x[2], max_x[2];				//dimensions of simulation region

	int max_list[2];						// contains max dimensions of domain

	vector<SPH_particle> particle_list;						//list of all the particles

	vector<vector<vector<SPH_particle*> > > search_grid;		//Outer 2 refers to the row and column of the grid respectively, inner vector is the list of pointers in each cell
};


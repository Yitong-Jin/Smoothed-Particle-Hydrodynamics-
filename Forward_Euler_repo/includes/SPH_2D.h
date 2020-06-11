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

	double x[2], v[2];				//position and velocity
	double rho, P;					//density and pressure
	double mass;
	double rho0;
	double dv[2];
	double drho;
	double visc;		// viscosity in Pa s
	double gamma;		// Stiffness
	bool boundary;
	int part_num;
	double rho_new;

	static SPH_main* main_data;		//link to SPH_main class so that it can be used in calc_index

	int list_num[2];				//index in neighbour finding array

	void calc_index(void);

	double get_B(void);

	double get_P(double& rho);

	double get_dwdr(double& r, double& h);

	void NS_acceleration_simultaneous(SPH_particle* other_particle, double* e_ij, double* v_ij, double& dist);

	void NS_acceleration(SPH_particle* other_particle, double* e_ij, double* v_ij, double dist);

	void NS_density_simultaneous(SPH_particle* other_particle, double* e_ij, double* v_ij, double& dist);

	void NS_density(SPH_particle* other_particle, double* e_ij, double* v_ij, double dist);

	double get_W(double& r, double& h);

	double repulsion_func(double& d, double& n);	// Defines a function that's repulsive at the boundaries and forceless elsewhere

	void apply_boundary_forces(string keyword);

	void update(string keyword);

	//void update(double dt);

	//~SPH_particle();
};


class SPH_main
{
public:
	SPH_main();

	void set_values(void);

	void initialise_grid(void);

	void place_points(double* min, double* max);

	void place_points_beach(double* min, double* max);

	void allocate_to_grid(void);			//allocates all the points to the search grid (assumes that index has been appropriately updated)

	void neighbour_iterate(SPH_particle* part);

	void smoothing();

	void smoothing_neighbour_itr(SPH_particle* part);

	void neighbour_iterate_stencil(SPH_particle* part);

	void set_adaptive_dt();

	double h_fac;							// set in set values - links h to dx^2 (1.3)

	double dx;								//particle initial spacing

	double g;

	double c_0;

	double dt;								// set in set_values
	double dt_CFL;
	double dt_F;
	double dt_A;

	double h;								//smoothing length

	double min_x[2], max_x[2];				//dimensions of simulation region

	int max_list[2];

	vector<SPH_particle> particle_list;						//list of all the particles

	vector<vector<vector<SPH_particle*> > > search_grid;		//Outer 2 are the grid, inner vector is the list of pointers in each cell
};
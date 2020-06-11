#include "SPH_2D.h"
#include "cmath"
#include "algorithm"
#include <string>
#define M_PI 3.14159265358979323846

// ==================================================================================
// ==================================================================================

//
// Initialising patricle class
//
SPH_particle::SPH_particle()
{
	/// Vector arrays
	
	///
	/// initiate x component of location
	///
	x[0][0] = 0.;	
	///
	/// initiate y component of location
	///
	x[1][0] = 0.;	
	///
	/// initiate x component of velocity
	///
	v[0][0] = 0.;						
	///
	/// initiate y component of velocity
	///
	v[1][0] = 0.;						
	rho[0] = 1000.;


	int step = 0;
	///
	/// Update gradient functions
	///
	dv[0] = 0.;
	dv[1] = 0.;
	drho = 0.;

	///
	/// Bulk properties
	///
	/// Initial density in kg/m^3
	///
	rho0 = 1000.;		
	///
	/// mass of particle
	///
	mass = pow(main_data->dx, 2) * rho0;	
	///
	/// viscosity in Pa s
	///
	visc = 0.001;		
	///
	/// Stiffness
	///
	gamma = 7.;	
	///
	/// density and pressure
	rho_new = 1000.;

	P = this->get_P(rho[step]);					
}

///
/// Initialise step for Improved Euler
///
void SPH_particle::set_step(int step_)
{
	this->step = step_;
}

///
/// Function to calculate B, where B is the reference density of water
///
double SPH_particle::get_B(void)
{
	return (this->rho0 * pow(main_data->c_0, 2.)) / this->gamma;
}

///
/// Function to get P, where P is pressure - calculated using the Tait Equation
///
double SPH_particle::get_P(double& rho)
{
	return get_B() * (pow((rho / this->rho0), gamma) - 1);
}

///
/// W is the smoothing kernel - has compact support and it is chosen to be normal in this SPH model
///
double SPH_particle::get_W(double& r, double& h) {

	double q = abs(r) / h;
	double C = 10. / (7. * M_PI * (h * h));
	if (q >= 0. && q <= 1.) {
		return C * (1. - 3. / 2. * pow(q, 2.) + 3. / 4. * pow(q, 3.));
	}
	else if (q > 1 && q <= 2.) {
		return C / 4. * pow((2. - q), 3.);
	}
	else {
		return 0.;
	}
}

///
/// Function to get the differential (dwdr) of the smoothing kernel
///
double SPH_particle::get_dwdr(double& r, double& h)
{
	double q = r / h;
	double C = 10. / (7. * M_PI * (h * h * h));
	if (q <= 1)
		return C * (-3. * q + 9. / 4. * q * q);
	else if (q <= 2.)
		return -3. * C / 4. * (2. - q) * (2. - q);
	else
		return 0.;
}

///
/// Function that calculates the acceleration of a particle based on summation using the smoothing kernel's gradient
/// NS_acceleration_simultaneous will be used in the neighbour_iterate_stencil and uses the fact that a pair of particles exert an equal and opposite force/ sum of the acceleration is anti-symmetric 
/// The time stepping constraint based on acceleration (dt_F) is also calculated in this function (adaptive time stepping used)
///
void SPH_particle::NS_acceleration_simultaneous(SPH_particle* other_particle, double* e_ij, double* v_ij, double& dist)
{
	///
	/// Initialise some values from the other particle's attributes
	///
	double op_rho = other_particle->rho[this->step];
	double op_mass = other_particle->mass;
	double op_visc = other_particle->visc;
	double op_P = other_particle->P;

	///
	/// Array to hold the acceleration in the x and y direction
	/// Initialise the inital values to be 0
	///
	double additive[2] = { 0, 0 };


	int other_step = (this->step + 1) % 2;

	///
	/// n = 0 is for the x-direction
	/// n = 1 is for the y-direction
	///

	for (int n = 0; n < 2; n++)
	{
		additive[n] += -op_mass * (this->P / pow(this->rho[other_step], 2) + op_P / pow(op_rho, 2)) * get_dwdr(dist, main_data->h) * e_ij[n];
		additive[n] += op_visc * op_mass * (1 / pow(this->rho[other_step], 2) + 1 / pow(op_rho, 2)) * get_dwdr(dist, main_data->h) * v_ij[n] / dist;

		this->dv[n] += additive[n];

		///
		/// if other_ particle is not in the boundary, calculate the acceleration
		///
		if (!other_particle->boundary)
			other_particle->dv[n] -= additive[n];
	}

	///
	/// Calculating temp_dt_F and only update dt_F if it is bigger than the temporary value
	///
	double tmp_F_dt = sqrt(main_data->h / sqrt((dv[0] * dv[0]) + (dv[1] * dv[1])));
	if (tmp_F_dt < main_data->dt_F) {
		main_data->dt_F = tmp_F_dt;
	}

}

///
/// Function that calculates the acceleration of a particle based on summation using the smoothing kernel's gradient
/// The time stepping constraint based on acceleration (dt_F) is also calculated in this function (adaptive time stepping used)
///
void SPH_particle::NS_acceleration(SPH_particle* other_particle, double* e_ij, double* v_ij, double dist)
{
	///
	/// Initialise some values from the other particle's attributes
	///
	double op_rho = other_particle->rho[this->step];
	double op_mass = other_particle->mass;
	double op_visc = other_particle->visc;
	double op_P = other_particle->P;

	///
	/// Accumulate the x-component of acceleration
	///
	int other_step = (this->step + 1) % 2;

	///
	/// n = 0 is for the x-direction
	/// n = 1 is for the y-direction
	///
	for (int n = 0; n < 2; n++)
	{
		this->dv[n] += -op_mass * (this->P / pow(this->rho[other_step], 2) + op_P / pow(op_rho, 2)) * get_dwdr(dist, main_data->h) * e_ij[n];
		this->dv[n] += op_visc * op_mass * (1 / pow(this->rho[other_step], 2) + 1 / pow(op_rho, 2)) * get_dwdr(dist, main_data->h) * v_ij[n] / dist;
	}

	///
	/// Calculating temp_dt_F and only update dt_F if it is bigger than the temporary value
	///
	double tmp_F_dt = sqrt(main_data->h / sqrt((dv[0] * dv[0]) + (dv[1] * dv[1])));
	if (tmp_F_dt < main_data->dt_F) {
		main_data->dt_F = tmp_F_dt;
	}
}

///
/// Function that calculates the density change of a particle over time based on summation using the smoothing kernel's gradient
/// NS_density_simultaneous will be used in the neighbour_iterate_stencil and uses the fact that the density interaction is also symmetric
/// The time stepping constraint based on density and acoustic speed (dt_A) is also calculated in this function (adaptive time stepping used)
///
void SPH_particle::NS_density_simultaneous(SPH_particle* other_particle, double* e_ij, double* v_ij, double& dist)
{
	int other_step = (this->step + 1) % 2;

	///
	/// Initialise the inital value to be 0
	///
	double additive = 0;
	///
	/// Calculate density addition
	///
	additive = other_particle->mass * get_dwdr(dist, main_data->h) * ((v_ij[0] * e_ij[0]) + (v_ij[1] * e_ij[1]));
	///
	/// Accumulate the density gradient
	///
	this->drho += additive;

	///
	/// This interaction is symmetric - so each increment of this drho can increment the other's drho
	///
	other_particle->drho += additive;

	///
	/// Calculating temp_dt_A and only update dt_A if it is bigger than the temporary value
	///
	double tmp_dt_A = main_data->h / (main_data->c_0 * sqrt(pow((rho[other_step] / rho0), gamma - 1)));

	if (tmp_dt_A < main_data->dt_A) {
		main_data->dt_A = tmp_dt_A;
	}
}

///
/// Function that calculates the density change of a particle over time based on summation using the smoothing kernel's gradient
/// The time stepping constraint based on density and acoustic speed (dt_A) is also calculated in this function (adaptive time stepping used)
///
void SPH_particle::NS_density(SPH_particle* other_particle, double* e_ij, double* v_ij, double dist)
{
	int other_step = (this->step + 1) % 2;
	///
	/// Accumulate the density gradient
	///
	this->drho += other_particle->mass * get_dwdr(dist, main_data->h) * ((v_ij[0] * e_ij[0]) + (v_ij[1] * e_ij[1]));

	///
	/// Calculating temp_dt_A and only update dt_A if it is bigger than the temporary value
	///
	double tmp_dt_A = main_data->h / (main_data->c_0 * sqrt(pow((rho[other_step] / rho0), gamma - 1)));


	if (tmp_dt_A < main_data->dt_A) {
		main_data->dt_A = tmp_dt_A;
	}
}

///
/// Function to calculate the repulsive force 
/// d is the distance between the particle and the boundary particle
/// Function contains an exponential based on d - for a constant value of n (2, empirically chosen), the repulsive force dramatically reduces as d increases
///
double SPH_particle::repulsion_func(double& d, double& n)
{
	return pow((0.5 * main_data->dx / d), n) * (1 - d / (0.5 * main_data->dx));
}

///
///Function to apply the repulsive forces onto the fluid particles
///
void SPH_particle::apply_boundary_forces(string keyword)
{
	///
	/// Initialise variables
	///

	///
	/// force of gravity
	///
	this->dv[1] -= main_data->g;	
	///
	/// the hydraulic head (height of fluid divided by 10, empirically chosen)
	///
	double h_dash = main_data->max_x[1] / 5;
	///
	/// pressure at a certain hydraulic head level is calculated by P = phg
	double P = this->rho0 * main_data->g * h_dash;
	///
	/// the power (n) for the repulsive force function
	///
	double n1 = 4;									
	///
	/// cut-off distance - if particle is within this cut-off distance from the boundary particle, then repulsive force will occur
	///
	double cut_off = 0.5 * main_data->dx;			
	
	///
	///In this model, right is positive (+) and upwards is positive (+)
	///
	///LEFT boundary
	///
	if (this->x[0][this->step] < cut_off) {
		///
		/// Add the repulsive forces to the acceleration for each of our dimensions
		///
		this->dv[0] += main_data->g * h_dash * repulsion_func(this->x[0][this->step], n1) / main_data->dx;
	}

	///
	/// RIGHT boundary
	///
	else if (this->x[0][this->step] > main_data->max_x[0] - 2 * main_data->h - cut_off) {
		double dist_to_edge = main_data->max_x[0] - 2 * main_data->h - this->x[0][this->step];
		///
		/// Add the repulsive forces to the acceleration for each of our dimensions
		///
		this->dv[0] -= main_data->g * h_dash * repulsion_func(dist_to_edge, n1) / main_data->dx;
	}

	///
	/// BOTTOM boundary
	///
	if (this->x[1][this->step] < cut_off) {
		///
		/// Add the repulsive forces to the acceleration for each of our dimensions
		///
		this->dv[1] += main_data->g * h_dash * repulsion_func(this->x[1][this->step], n1) / main_data->dx;
	}


	// TOP boundary
	else if (this->x[1][this->step] > main_data->max_list[1] - 2 * main_data->h - cut_off) {
		double dist_to_edge = main_data->max_x[1] - 2 * main_data->h - this->x[1][this->step];
		// Add the repulsive forces to the acceleration for each of our dimensions
		this->dv[1] -= main_data->g * h_dash * repulsion_func(dist_to_edge, n1) / main_data->dx;
	}

	if (keyword == "beach")
	{
		///
		/// BOTTOM-beach
		///
		if ((this->x[0][this->step] >= (main_data->max_x[0] - 2 * main_data->h) * (7.0 / 8.0)) && ((this->x[1][this->step] - ((main_data->max_x[0] - 2 * main_data->h) * (1.0 / 8.0))) < cut_off))
		{
			///
			/// Calculate distance to edge of domain
			///
			double dist_to_edge = this->x[1][this->step] - ((main_data->max_x[0] - 2 * main_data->h) * (1.0 / 8.0));
			this->dv[1] += main_data->g * h_dash * repulsion_func(dist_to_edge, n1) / main_data->dx;
		}
		///
		/// apply slope
		///
		if ((this->x[0][this->step] >= (main_data->max_x[0] - 2 * main_data->h) * (3.0 / 4.0)) && (this->x[0][this->step] < (main_data->max_x[0] - 2 * main_data->h) * (7.0 / 8.0)))
		{
			///
			/// define slope points
			///
			double slope_x0 = this->x[0][this->step] - (main_data->max_x[0] - 2 * main_data->h) * (3.0 / 4.0);
			double dist_to_edge = this->x[1][this->step] - slope_x0;

			if (((this->x[1][this->step] - slope_x0) < cut_off) && (dist_to_edge < cut_off))
			{
				///
				/// Add the repulsive forces to the acceleration for each of our dimensions
				///
				this->dv[0] -= main_data->g * h_dash * repulsion_func(dist_to_edge, n1) / main_data->dx;
				this->dv[1] += main_data->g * h_dash * repulsion_func(dist_to_edge, n1) / main_data->dx;
			}
		}
	}
}

///
/// Function to update the properties of the particles after every iteration (deprecated, see Forward Euler Repo)
///
void SPH_particle::update(string keyword)
{
	this->step = 0;

	///
	/// Update the time-step to use an appropriate one for the current iteration
	///
	main_data->set_adaptive_dt();

	///
	/// Apply repulsive force only to non-boundary particles
	///
	if (!this->boundary)
		this->apply_boundary_forces(keyword);

	///
	/// Update the position and the velocities for non-boundary particles
	///
	double dt = this->main_data->dt;
	if (!this->boundary) {
		for (int n = 0; n < 2; n++)
		{
			this->x[n][this->step] += dt * this->v[n][this->step];
			this->v[n][this->step] += dt * this->dv[n];
		}
	}

	///
	/// Update the density
	///
	this->rho[this->step] += dt * this->drho;
	
	///
	/// Update the pressure
	///
	this->P = get_P(this->rho[this->step]);

	///
	/// Reset gradients (in boundaries as well to prevent overflow)
	///
	for (int n = 0; n < 2; n++)
		this->dv[n] = 0.0;
	this->drho = 0.;

	///
	/// Update the row and col index of the particle's sub-grid
	///
	this->calc_index();
}

///
/// Function to update particle for both Forward and Improved Euler
///
void SPH_particle::update_pc(int step, string keyword)
{
	this->step = step;
	main_data->set_adaptive_dt();
	if (!this->boundary)
		this->apply_boundary_forces(keyword);

	double dt = this->main_data->dt;

	if (step == 1)
	{
		///
		/// if not a boundary particle
		///
		if (!this->boundary) {						
			for (int n = 0; n < 2; n++)
			{
				this->x[n][step] = this->x[n][0] + 0.5 * dt * this->v[n][0];
				this->v[n][step] = this->v[n][0] + 0.5 * dt * this->dv[n];
			}
		}
		
		///
		/// Update the density
		///
		this->rho[step] = this->rho[0] + 0.5 * dt * this->drho;

	}
	else
	{
		if (!this->boundary)
		{
			for (int n = 0; n < 2; n++)
			{

				this->x[n][0] = x[n][0] + main_data->dt * this->v[n][1];
				this->v[n][0] += main_data->dt * this->dv[n];
			}
		}

		this->rho[0] = this->rho[0] + dt * this->drho;
	}

	///
	/// Update the pressure
	///
	this->P = get_P(this->rho[step]);

	///
	/// Reset gradients (in boundaries as well to prevent overflow)
	///
	for (int n = 0; n < 2; n++)
		this->dv[n] = 0.0;
	this->drho = 0.;

	///
	/// Update the row and col index of the particle's sub-grid
	///
	this->calc_index();
}

///
/// Link pointer instance of SPH_main to SPH_particle's main data
///
SPH_main* SPH_particle::main_data;

///
/// Function to calculate the row and col index of the particle's sub-grid
///
void SPH_particle::calc_index(void)
{
	for (int i = 0; i < 2; i++)
		list_num[i] = int((x[i][this->step] - main_data->min_x[i]) / (2.0 * main_data->h));
}

///
/// Function that smooths densities of all the particles
///
void SPH_main::smoothing()
{
	///
	/// Smooths all particles
	///
	for (int i = 0; i < particle_list.size(); i++)
	{
		smoothing_neighbour_itr(&particle_list[i]);
	}

	///
	/// Updates the density and pressure
	///
	for (int i = 0; i < particle_list.size(); i++)
	{
		particle_list[i].rho[particle_list[i].step] = particle_list[i].rho_new;
		particle_list[i].P = particle_list[i].get_P(particle_list[i].rho[particle_list[i].step]);
	}
}

///
/// Function to smooth the density for a particular particle
/// Requires iterating between the particle's neighbours
///
void SPH_main::smoothing_neighbour_itr(SPH_particle* part)
{
	///
	/// Initialise neighbouring particle
	///
	SPH_particle* other_part;
	///
	/// Euclidean distance between particles
	///
	double dist;
	///
	/// array to hold the x and y distances between particles
	///
	double dn[2];			
	double W_sum = 0.;
	double W_sum_rho = 0.;
	///
	/// neighbour count - if particle has no neighbours within 2h, then no need to smooth density
	///
	int nb_neighbour = 0;   
	part->rho_new = part->rho[part->step];

	//int other_step = (part->step + 1) % 2;

	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0])
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1]) {
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++) {
						other_part = search_grid[i][j][cnt];

						///
						///Calculates the distance between potential neighbours
						///
						for (int n = 0; n < 2; n++)
							dn[n] = part->x[n][part->step] - other_part->x[n][other_part->step];
						dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

						///
						/// only particle within 2h
						///
						if (dist < 2. * h)					
						{
							double w = part->get_W(dist, this->h);
							W_sum += w;
							W_sum_rho += w / other_part->rho[other_part->step];
							nb_neighbour++;
						}
					}
				}

	///
	/// Equation to smooth the density (if there are neighbours)
	///
	if (nb_neighbour > 0) { part->rho_new = W_sum / W_sum_rho; }
}

///
/// Function to update the time-step according to the properties of the current time
///
void SPH_main::set_adaptive_dt()
{
	///
	/// the time-step will be equal to the smallest of the three constraints
	///
	this->dt = 0.1 * min({ dt_CFL, dt_F, dt_A });
}

///
/// Initialise so that SPH_particle main_data points to SPH_main
///
SPH_main::SPH_main()
{
	SPH_particle::main_data = this;
}

///
/// Set values for the entire space
///
void SPH_main::set_values(void)
{
	min_x[0] = 0.0;
	min_x[1] = 0.0;

	///
	/// Space will be 20m wide
	max_x[0] = 20.0; 
	///
	/// Space will be 10m in height
	///
	max_x[1] = 10.0; 

	dx = 0.2;

	c_0 = 20;

	g = 9.81;

	h_fac = 1.3;

	h = dx * h_fac;

	dt = 0.1 * h / c_0;
	dt_CFL = 10e8;
	dt_F = 10e8;
	dt_A = 10e8;
}

///
/// Function to initialise a grid onto the space
///
void SPH_main::initialise_grid(void)
{
	///
	/// calculating the true size of the grid - add buffer for virtual wall particles
	///
	for (int i = 0; i < 2; i++)
	{
		min_x[i] -= 2.0 * h;
		max_x[i] += 2.0 * h;

		max_list[i] = int((max_x[i] - min_x[i]) / (2.0 * h) + 1.0);
	}

	///
	/// initialise the number of rows
	///
	search_grid.resize(max_list[0]);
	for (int i = 0; i < max_list[0]; i++)
		///
		/// initialise the number of columns
		///
		search_grid[i].resize(max_list[1]);
}

///
/// Function to place both the fluid and boundary particles into the grid
///
void SPH_main::place_points(double* min, double* max)
{
	///
	/// initialise an array to hold the starting values for the boundaries
	///
	double x[2] = { min[0], min[1] };
	SPH_particle particle;
	particle.step = 0;

	///
	/// place bottom boundary and add particles to particle list, plus calculate subgrid index
	///
	x[0] = 0;
	while (x[0] <= (max[0] - 2 * h))
	{
		x[1] = min[1];
		while (x[1] <= 0)
		{
			for (int i = 0; i < 2; i++)
				particle.x[i][particle.step] = x[i];

			particle.calc_index();
			particle.boundary = true;
			particle_list.push_back(particle);

			x[1] += dx;
		}
		x[0] += dx;
	}

	///
	/// place left boundary and add particles to particle list, plus calculate subgrid index
	///
	x[0] = min[0];
	x[1] = min[1];
	while (x[0] <= 0)
	{
		x[1] = min[1];
		while (x[1] <= max[1])
		{
			for (int i = 0; i < 2; i++)
				particle.x[i][particle.step] = x[i];

			particle.calc_index();
			particle.boundary = true;
			particle_list.push_back(particle);

			x[1] += dx;
		}
		x[0] += dx;
	}


	///
	/// place right boundary and add particles to particle list, plus calculate subgrid index
	///
	x[0] = max[0] - 2.0 * h;
	x[1] = min[0];
	while (x[0] <= max[0])
	{
		x[1] = min[0];
		while (x[1] <= max[1])
		{
			for (int i = 0; i < 2; i++)
				particle.x[i][particle.step] = x[i];

			particle.calc_index();
			particle.boundary = true;
			particle_list.push_back(particle);

			x[1] += dx;
		}
		x[0] += dx;
	}

	///
	/// place top boundary and add particles to particle list, plus calculate subgrid index
	///
	x[0] = 0;
	x[1] = max[1] - 1.7 * h;
	while (x[0] <= (max[0] - 2.0 * h))
	{
		x[1] = max[1] - 1.7 * h;
		while (x[1] <= max[1])
		{
			for (int i = 0; i < 2; i++)
				particle.x[i][particle.step] = x[i];

			particle.calc_index();
			particle.boundary = true;
			particle_list.push_back(particle);

			x[1] += dx;
		}
		x[0] += dx;
	}

	///
	/// place fluid particles
	///
	///
	/// create arrays to hold the height and width of the water body
	/// array that holds the range x-values for the leftmost points of the lake
	///
	double init_x[2] = { dx / 2, dx / 2 };
	///
	/// array that holds the range of y-values for the leftmost points of the lake
	///
	double init_y[2] = { dx / 2, 2.0 + dx / 2 };
	///
	/// array that holds the range x-values for the rightmost points of the lake
	///
	double set_x[2] = { max[0] - 2.0 * h, 3 };
	///
	/// array that holds the range y-values for the rightmost points of the lake
	///
	double set_y[2] = { 2, 5 };

	///
	/// loop over the arrays and increment using particle spacing dx to place the fluid particles into the grid
	///
	for (int i = 0; i < 2; i++)
	{
		x[0] = init_x[i];
		x[1] = init_y[i];
		while (x[0] < set_x[i])
		{
			x[1] = init_y[i];
			while (x[1] < set_y[i])
			{
				for (int i = 0; i < 2; i++)
					particle.x[i][particle.step] = x[i];

				particle.boundary = false;
				particle.calc_index();

				particle_list.push_back(particle);

				x[1] += dx;
			}
			x[0] += dx;
		}
	}

}

///
/// Function to place particles in the beach domain
///
void SPH_main::place_points_beach(double* min, double* max)
{
	double x[2] = { min[0], min[1] };
	SPH_particle particle;
	particle.step = 0;
	double temp_pp;

	///
	/// place bottom boundary
	/// bottom-left
	///
	x[0] = 0;
	while (x[0] <= (max[0] - 2 * h) * 3.0 / 4.0)
	{
		x[1] = min[1];
		while (x[1] <= 0)
		{
			for (int i = 0; i < 2; i++)
				particle.x[i][particle.step] = x[i];

			particle.calc_index();
			particle.boundary = true;
			particle_list.push_back(particle);

			x[1] += dx;
		}
		x[0] += dx;
	}

	///
	/// bottom-right-beach
	///
	x[0] = (max[0] - 2 * h) * (7.0 / 8.0);
	while (x[0] <= (max[0] - 2 * h))
	{
		x[1] = (max[0] - 2 * h) * (1.0 / 8.0) - 2 * h;
		while (x[1] <= (max[0] - 2 * h) * (1.0 / 8.0))
		{
			for (int i = 0; i < 2; i++)
				particle.x[i][particle.step] = x[i];

			particle.calc_index();
			particle.boundary = true;
			particle_list.push_back(particle);

			x[1] += dx;
		}
		x[0] += dx;
	}

	///
	/// bottom-slope
	///
	x[0] = (max[0] - 2 * h) * (3.0 / 4.0);
	temp_pp = min[1];
	while (x[0] <= ((max[0] - 2 * h) * (7.0 / 8.0)))
	{
		x[1] = temp_pp;
		while (x[1] <= (temp_pp + 2 * h))
		{
			for (int i = 0; i < 2; i++)
				particle.x[i][particle.step] = x[i];

			particle.calc_index();
			particle.boundary = true;
			particle_list.push_back(particle);

			x[1] += dx;
		}
		temp_pp += dx;
		x[0] += dx;
	}

	///
	/// place left boundary
	///
	x[0] = min[0];
	x[1] = min[0];
	while (x[0] <= 0)
	{
		x[1] = min[0];
		while (x[1] <= max[1])
		{
			for (int i = 0; i < 2; i++)
				particle.x[i][particle.step] = x[i];

			particle.calc_index();
			particle.boundary = true;
			particle_list.push_back(particle);

			x[1] += dx;
		}
		x[0] += dx;
	}

	///
	/// place right boundary
	///
	x[0] = max[0] - 2.0 * h;
	x[1] = (max[0] - 2 * h) * (1.0 / 8.0) - 2 * h;
	while (x[0] <= max[0])
	{
		x[1] = (max[0] - 2 * h) * (1.0 / 8.0) - 2 * h;
		while (x[1] <= max[1])
		{
			for (int i = 0; i < 2; i++)
				particle.x[i][particle.step] = x[i];

			particle.calc_index();
			particle.boundary = true;
			particle_list.push_back(particle);

			x[1] += dx;
		}
		x[0] += dx;
	}

	///
	/// place top boundary
	///
	x[0] = 0;
	x[1] = max[1] - 1.7 * h;
	while (x[0] <= (max[0] - 2.0 * h))
	{
		x[1] = max[1] - 1.7 * h;
		while (x[1] <= max[1])
		{
			for (int i = 0; i < 2; i++)
				particle.x[i][particle.step] = x[i];

			particle.calc_index();
			particle.boundary = true;
			particle_list.push_back(particle);

			x[1] += dx;
		}
		x[0] += dx;
	}

	///
	/// place fluid particle
	///
	double init_x[2] = { dx / 2, dx / 2 };
	double init_y[2] = { dx / 2, 2.0 + dx / 2 };
	double set_x[2] = { (max[0] - 2.0 * h) * 3.0 / 4.0, 3 };
	double set_y[2] = { 2, 5 };

	for (int i = 0; i < 2; i++)
	{
		x[0] = init_x[i];
		x[1] = init_y[i];
		while (x[0] < set_x[i])
		{
			x[1] = init_y[i];
			while (x[1] < set_y[i])
			{
				for (int i = 0; i < 2; i++)
					particle.x[i][particle.step] = x[i];

				particle.boundary = false;
				particle.calc_index();

				particle_list.push_back(particle);

				x[1] += dx;
			}
			x[0] += dx;
		}
	}
	///
	/// slope fluid particle
	///
	x[0] = (max[0] - 2 * h) * (3.0 / 4.0);// +dx / 2;
	temp_pp = dx / 2;
	while (x[0] <= ((max[0] - 2 * h) * (3.0 / 4.0) + set_y[0]))
	{
		x[1] = temp_pp;
		while (x[1] <= set_y[0])
		{
			for (int i = 0; i < 2; i++)
				particle.x[i][particle.step] = x[i];

			particle.calc_index();
			particle.boundary = false;
			particle_list.push_back(particle);

			x[1] += dx;
		}
		temp_pp += dx;
		x[0] += dx;
	}
}

///
/// Function that clears and updates the grid with particles at the current time-step
/// Thus, the function needs to be called after each iteration so that all the particles have their positions updated
///
void SPH_main::allocate_to_grid(void)
{
	///
	/// clear the old grid
	///
	for (int i = 0; i < max_list[0]; i++)
		for (int j = 0; j < max_list[1]; j++)
			search_grid[i][j].clear();

	for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
	{
		///
		/// update the particle number according to which subgrid its in
		///
		particle_list[cnt].part_num = search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].size();
		///
		/// add the particle to its corresponding subgrid
		///
		search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
	}
}

///
/// Function to calculate the acceleration and density by iterating over all particles within 2h of the chosen particle 
/// Iteration occurs for all 8 neighbouring sub-grids in a step-by-step manner
/// Also updates dt_CFL
///
void SPH_main::neighbour_iterate(SPH_particle* part)
{
	///
	/// Initialise other particle
	///
	SPH_particle* other_part;
	///
	/// distance between particles
	///
	double dist;			
	///
	/// array to hold the x and y distances between particles
	///
	double dn[2];			
	double tmp_CFL_dt;
	double e_ij[2], v_ij[2];

	double f_rep = 0;

	int other_step = (part->step + 1) % 2;

	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0]) {
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1])
				{
					///
					/// calculate dist_to_edge x and y for particle
					///
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
					{
						///
						/// re-direct the address
						///
						other_part = search_grid[i][j][cnt];	
						///
						///stops particle interacting with itself
						///

						if (part != other_part)					
						{
							///
							/// Calculates the distance between potential neighbours
							///
							for (int n = 0; n < 2; n++)
								dn[n] = part->x[n][other_step] - other_part->x[n][other_step];

							dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

							///
							///only particle within 2h
							///
							if (dist < 2. * h && dist > 0.)
							{
								///
								///all the interactions between the particles
								///
								for (int n = 0; n < 2; n++)
								{
									///
									/// difference in velocity
									///
									v_ij[n] = part->v[n][other_step] - other_part->v[n][other_step];
									///
									/// difference in direction
									///
									e_ij[n] = (part->x[n][other_step] - other_part->x[n][other_step]) / dist;
								}
								
								///
								/// calculate the Courant number to check timestep is still valid
								///
								tmp_CFL_dt = h / sqrt((v_ij[0] * v_ij[0]) + (v_ij[1] * v_ij[1]));
								if (tmp_CFL_dt < this->dt_CFL) {
									this->dt_CFL = tmp_CFL_dt;
								}

								///
								/// calculate particle's acceleration
								///
								part->NS_acceleration(other_part, e_ij, v_ij, dist);
								///
								/// calculate particle's density
								///
								part->NS_density(other_part, e_ij, v_ij, dist);
							}
						}
					}
				}
		}
}

///
/// Function to calculate the acceleration and density by iterating over particles within 2h of the chosen particle
/// Unlike the function above (non-stencil), this function iterates only for 4 neighbouring sub-grids (in an L-shape)
/// This can be used because the interations are symmetric, so interations for the 4 neighbouring sub-grids will be the same as the other four
/// More efficient as halves number of calculations
/// Also updates dt_CFL
///
void SPH_main::neighbour_iterate_stencil(SPH_particle* part)
{
	///
	/// Initialise other particle
	///
	SPH_particle* other_part;
	///
	/// distance between particles
	///
	double dist;			
	///
	/// array to hold the x and y distances between particles
	///
	double dn[2];			

	///
	/// Array holding the row index of the 4 neighbouring sub-grids and its own grid
	///
	signed int row_nei[5] = { 0, 1, -1, 0, 1 };
	///
	// Array holding the column index of the 4 neighbouring sub-grids and its own grid
	///
	signed int col_nei[5] = { 0, 0,  1, 1, 1 };

	double tmp_CFL_dt;

	int other_step = (part->step + 1) % 2;
	double e_ij[2], v_ij[2];

	for (int k = 0; k < 5; k++)
	{
		int i = part->list_num[0] + row_nei[k];
		int j = part->list_num[1] + col_nei[k];
		if (i >= 0 && i < max_list[0] && j >= 0 && j < max_list[1])
		{
			if (k != 0)
			{
				for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
				{
					other_part = search_grid[i][j][cnt];

					///
					/// Calculates the distance between potential neighbours
					///
					for (int n = 0; n < 2; n++)
						dn[n] = part->x[n][other_step] - other_part->x[n][other_step];

					dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

					///
					/// only particle within 2h
					///
					if (dist < 2. * h && dist > 0.)
					{
						///
						/// All the interactions between the particles
						///
						for (int n = 0; n < 2; n++)
						{
							v_ij[n] = part->v[n][other_step] - other_part->v[n][other_step];
							e_ij[n] = (part->x[n][other_step] - other_part->x[n][other_step]) / dist;
						}

						part->NS_acceleration_simultaneous(other_part, e_ij, v_ij, dist);
						part->NS_density_simultaneous(other_part, e_ij, v_ij, dist);
					}
				}
			}
			else
			{
				//cout << "search_grid[i][j].size(): " << search_grid[i][j].size() << 
				for (unsigned int cnt = part->part_num + 1; cnt < search_grid[i][j].size(); cnt++)
				{
					other_part = search_grid[i][j][cnt];

					///
					/// Calculates the distance between potential neighbours
					///
					for (int n = 0; n < 2; n++)
						dn[n] = part->x[n][other_step] - other_part->x[n][other_step];

					dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

					///
					/// only particle within 2h
					///
					if (dist < 2. * h && dist > 0.)					
					{
						///
						// All the interactions between the particles
						///
						for (int n = 0; n < 2; n++)
						{
							v_ij[n] = part->v[n][other_step] - other_part->v[n][other_step];
							e_ij[n] = (part->x[n][other_step] - other_part->x[n][other_step]) / dist;
						}

						///
						/// Calculate the Courant number to check timestep is still valid
						///
						tmp_CFL_dt = h / sqrt((v_ij[0] * v_ij[0]) + (v_ij[1] * v_ij[1]));
						if (tmp_CFL_dt < this->dt_CFL) {
							this->dt_CFL = tmp_CFL_dt;
						}
						
						///
						/// calculate particle's acceleration
						///
						part->NS_acceleration_simultaneous(other_part, e_ij, v_ij, dist);
						///
						/// calculate particle's density
						///
						part->NS_density_simultaneous(other_part, e_ij, v_ij, dist);
					}

				}
			}
		}
	}
}

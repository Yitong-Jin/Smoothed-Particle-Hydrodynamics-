#include "../includes/SPH_2D.h"
#include "cmath"
#include "algorithm"
#define M_PI 3.14159265358979323846

// ==================================================================================
// ==================================================================================

SPH_particle::SPH_particle()
{
	// Vector arrays
	x[0] = 0.;						// initiate x component of location
	x[1] = 0.;						// initiate y component of location
	v[0] = 0.;						// initiate x component of velocity
	v[1] = 0.;						// initiate y component of velocity
	rho = 1000.;


	// Update gradient functions
	dv[0] = 0.;
	dv[1] = 0.;
	drho = 0.;


	// Bulk properties
	rho0 = 1000.;							// Initial density in kg/m^3
	mass = pow(main_data->dx, 2) * rho0;	// mass of particle
	visc = 0.001;							// viscosity in Pa s
	gamma = 7.;								// Stiffness
	rho_new = 1000.;

	P = this->get_P(rho);					//density and pressure
}

double SPH_particle::get_B(void)
{
	return (this->rho0 * pow(main_data->c_0, 2.)) / this->gamma;
}

double SPH_particle::get_P(double& rho)
{
	return get_B() * (pow((rho / this->rho0), gamma) - 1);
}

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

void SPH_particle::NS_acceleration_simultaneous(SPH_particle* other_particle, double* e_ij, double* v_ij, double& dist)
{
	// Initialise some values from the other particle's attributes
	double op_rho = other_particle->rho;
	double op_mass = other_particle->mass;
	double op_visc = other_particle->visc;
	double op_P = other_particle->P;

	double additive[2] = { 0, 0 };

	for (int n = 0; n < 2; n++)
	{
		additive[n] += -op_mass * (this->P / pow(this->rho, 2) + op_P / pow(op_rho, 2)) * get_dwdr(dist, main_data->h) * e_ij[n];
		additive[n] += op_visc * op_mass * (1 / pow(this->rho, 2) + 1 / pow(op_rho, 2)) * get_dwdr(dist, main_data->h) * v_ij[n] / dist;

		this->dv[n] += additive[n];
		if (!other_particle->boundary)
			other_particle->dv[n] -= additive[n];
	}

	double tmp_F_dt = sqrt(main_data->h / sqrt((dv[0] * dv[0]) + (dv[1] * dv[1])));
	if (tmp_F_dt < main_data->dt_F) {
		main_data->dt_F = tmp_F_dt;
	}

}

void SPH_particle::NS_acceleration(SPH_particle* other_particle, double* e_ij, double* v_ij, double dist)
{
	// Initialise some values from the other particle's attributes
	double op_rho = other_particle->rho;
	double op_mass = other_particle->mass;
	double op_visc = other_particle->visc;
	double op_P = other_particle->P;
	//cout << ex_ij << " " << ey_ij << " " << vx_ij << " " << vy_ij << endl;
	// Accumulate the x-component of acceleration

	for (int n = 0; n < 2; n++)
	{
		this->dv[n] += -op_mass * (this->P / pow(this->rho, 2) + op_P / pow(op_rho, 2)) * get_dwdr(dist, main_data->h) * e_ij[n];
		this->dv[n] += op_visc * op_mass * (1 / pow(this->rho, 2) + 1 / pow(op_rho, 2)) * get_dwdr(dist, main_data->h) * v_ij[n] / dist;
	}

	double tmp_F_dt = sqrt(main_data->h / sqrt((dv[0] * dv[0]) + (dv[1] * dv[1])));
	if (tmp_F_dt < main_data->dt_F) {
		main_data->dt_F = tmp_F_dt;
	}
}

void SPH_particle::NS_density_simultaneous(SPH_particle* other_particle, double* e_ij, double* v_ij, double& dist)
{
	double additive = 0;
	// Calculate density addition
	additive = other_particle->mass * get_dwdr(dist, main_data->h) * ((v_ij[0] * e_ij[0]) + (v_ij[1] * e_ij[1]));
	// Accumulate the density gradient
	this->drho += additive;
	// This interaction is symmetric - so each increment of this drho can increment the other's drho

	other_particle->drho += additive;
	double tmp_dt_A = main_data->h / (main_data->c_0 * sqrt(pow((rho / rho0), gamma - 1)));
	if (tmp_dt_A < main_data->dt_A) {
		main_data->dt_A = tmp_dt_A;
	}
}

void SPH_particle::NS_density(SPH_particle* other_particle, double* e_ij, double* v_ij, double dist)
{
	// Accumulate the density gradient
	this->drho += other_particle->mass * get_dwdr(dist, main_data->h) * ((v_ij[0] * e_ij[0]) + (v_ij[1] * e_ij[1]));

	double tmp_dt_A = main_data->h / (main_data->c_0 * sqrt(pow((rho / rho0), gamma - 1)));
	if (tmp_dt_A < main_data->dt_A) {
		main_data->dt_A = tmp_dt_A;
	}
}

double SPH_particle::repulsion_func(double& d, double& n)
{
	return pow((0.5 * main_data->dx / d), n) * (1 - d / (0.5 * main_data->dx));
}

void SPH_particle::apply_boundary_forces(string keyword)
{

	this->dv[1] -= main_data->g;
	double h_dash = main_data->max_x[1] / 10;			// approx. hydraulic head
	double P = this->rho0 * main_data->g * h_dash;
	double n1 = 2;
	double cut_off = 0.5 * main_data->dx;

	// POSITIVE RIGHT  +

	// POSITIVE UPWARD +
	// if either of our particles are boundaries, and the distance between the particles is less than cut-off distance
	// LEFT
	if (this->x[0] < cut_off) {
		// Add the repulsive forces to the acceleration for each of our dimensions
		this->dv[0] += main_data->g * h_dash * repulsion_func(this->x[0], n1) / main_data->dx;
	}
	// if either of our particles are boundaries, and the distance between the particles is less than cut-off distance
	// RIGHT
	else if (this->x[0] > main_data->max_x[0] - 2 * main_data->h - cut_off) {
		double dist_to_edge = main_data->max_x[0] - 2 * main_data->h - this->x[0];
		// Add the repulsive forces to the acceleration for each of our dimensions
		this->dv[0] -= main_data->g * h_dash * repulsion_func(dist_to_edge, n1) / main_data->dx;
	}
	// BOTTOM
	if (this->x[1] < cut_off) {
		// Add the repulsive forces to the acceleration for each of our dimensions
		this->dv[1] += main_data->g * h_dash * repulsion_func(this->x[1], n1) / main_data->dx;
	}
	// TOP
	else if (this->x[1] > main_data->max_list[1] - 2 * main_data->h - cut_off) {
		double dist_to_edge = main_data->max_x[1] - 2 * main_data->h - this->x[1];
		// Add the repulsive forces to the acceleration for each of our dimensions
		this->dv[1] -= main_data->g * h_dash * repulsion_func(dist_to_edge, n1) / main_data->dx;
	}

	if (keyword == "beach")
	{
		// BOTTOM-beach
		if ((this->x[0] >= (main_data->max_x[0] - 2 * main_data->h) * (7.0 / 8.0)) && ((this->x[1] - ((main_data->max_x[0] - 2 * main_data->h) * (1.0 / 8.0))) < cut_off))
		{
			double dist_to_edge = this->x[1] - ((main_data->max_x[0] - 2 * main_data->h) * (1.0 / 8.0));
			this->dv[1] += main_data->g * h_dash * repulsion_func(dist_to_edge, n1) / main_data->dx;
		}
		// slope
		if ((this->x[0] >= (main_data->max_x[0] - 2 * main_data->h) * (3.0 / 4.0)) && (this->x[0] < (main_data->max_x[0] - 2 * main_data->h) * (7.0 / 8.0)))
		{
			double slope_x0 = this->x[0] - (main_data->max_x[0] - 2 * main_data->h) * (3.0 / 4.0);
			double dist_to_edge = this->x[1] - slope_x0;

			if (((this->x[1] - slope_x0) < cut_off) && (dist_to_edge < cut_off))
			{
				this->dv[0] -= main_data->g * h_dash * repulsion_func(dist_to_edge, n1) / main_data->dx;
				this->dv[1] += main_data->g * h_dash * repulsion_func(dist_to_edge, n1) / main_data->dx;
			}
		}
	}
}

void SPH_particle::update(string keyword)
{
	main_data->set_adaptive_dt();

	if (!this->boundary)
		this->apply_boundary_forces(keyword);

	double dt = this->main_data->dt;
	if (!this->boundary) {						// if not a boundary particle
		for (int n = 0; n < 2; n++)
		{
			this->x[n] += dt * this->v[n];
			this->v[n] += dt * this->dv[n];
		}
	}
	this->rho += dt * this->drho;
	this->P = get_P(this->rho);
	// Reset gradients (in boundaries as well to prevent overflow)
	for (int n = 0; n < 2; n++)
		this->dv[n] = 0.0;
	this->drho = 0.;
	this->calc_index();
}

SPH_main* SPH_particle::main_data;

void SPH_particle::calc_index(void)
{
	for (int i = 0; i < 2; i++)
		list_num[i] = int((x[i] - main_data->min_x[i]) / (2.0 * main_data->h));
}

void SPH_main::smoothing()
{
	for (int i = 0; i < particle_list.size(); i++)
	{
		smoothing_neighbour_itr(&particle_list[i]);
	}

	for (int i = 0; i < particle_list.size(); i++)
	{
		particle_list[i].rho = particle_list[i].rho_new;
		particle_list[i].P = particle_list[i].get_P(particle_list[i].rho);
	}
}

void SPH_main::smoothing_neighbour_itr(SPH_particle* part)
{
	SPH_particle* other_part;
	double dist;			//distance between particles
	double dn[2];			//vector from 1st to 2nd particle
	double W_sum = 0.;
	double W_sum_rho = 0.;
	int nb_neighbour = 0;
	part->rho_new = part->rho;

	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0])
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1]) {
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++) {
						other_part = search_grid[i][j][cnt];

						//Calculates the distance between potential neighbours
						for (int n = 0; n < 2; n++)
							dn[n] = part->x[n] - other_part->x[n];
						dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

						if (dist < 2. * h)					//only particle within 2h
						{
							//TODO: all the interactions between the particles
							double w = part->get_W(dist, this->h);
							W_sum += w;
							W_sum_rho += w / other_part->rho;
							nb_neighbour++;
						}
					}
				}

	if (nb_neighbour > 0) { part->rho_new = W_sum / W_sum_rho; }
}

void SPH_main::set_adaptive_dt()
{
	double min_dt = 0.15 * min({ dt_CFL, dt_F, dt_A });
	if (min_dt < this->dt) { this->dt = min_dt; }
}

SPH_main::SPH_main()
{
	SPH_particle::main_data = this;
}

void SPH_main::set_values(void)
{
	min_x[0] = 0.0;
	min_x[1] = 0.0;

	max_x[0] = 20.0;
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

void SPH_main::initialise_grid(void)
{
	for (int i = 0; i < 2; i++)
	{
		min_x[i] -= 2.0 * h;
		max_x[i] += 2.0 * h;												//add buffer for virtual wall particles

		max_list[i] = int((max_x[i] - min_x[i]) / (2.0 * h) + 1.0);
	}

	search_grid.resize(max_list[0]);
	for (int i = 0; i < max_list[0]; i++)
		search_grid[i].resize(max_list[1]);
}

void SPH_main::place_points(double* min, double* max)
{
	double x[2] = { min[0], min[1] };
	SPH_particle particle;
	// place boundary particles
	// array that holds the range x-values for the leftmost points of the boundary
	double init_x_b[4] = { 0, min[0], max[0] - 2.0 * h, 0. };
	// array that holds the range of y-values for the leftmost points of the boundary
	double init_y_b[4] = { min[1], min[1], min[0], max[1] - 1.7 * h };
	// array that holds the range x-values for the rightmost points of the boundary
	double set_x_b[4] = { max[0] - 2 * h, 0, max[0], max[0] - 2.0 * h };
	// array that holds the range y-values for the rightmost points of the boundary
	double set_y_b[4] = { 0,  max[1],  max[1], max[1] };
	for (int i = 0; i < 4; i++)
	{
		x[0] = init_x_b[i];
		x[1] = init_y_b[i];
		while (x[0] < set_x_b[i])
		{
			x[1] = init_y_b[i];
			while (x[1] < set_y_b[i])
			{
				for (int i = 0; i < 2; i++)
					particle.x[i] = x[i];
				particle.boundary = true;
				particle.calc_index();
				particle_list.push_back(particle);
				x[1] += dx;
			}
			x[0] += dx;
		}
	}
	// place fluid particle
	// array that holds the range x-values for the leftmost points of the lake
	double init_x[2] = { dx / 2, dx / 2 };
	// array that holds the range of y-values for the leftmost points of the lake
	double init_y[2] = { dx / 2, 2.0 + dx / 2 };
	// array that holds the range x-values for the rightmost points of the boundary
	double set_x[2] = { max[0] - 2.0 * h, 3 };
	// array that holds the range y-values for the rightmost points of the lake
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
					particle.x[i] = x[i];
				particle.boundary = false;
				particle.calc_index();
				particle_list.push_back(particle);
				x[1] += dx;
			}
			x[0] += dx;
		}
	}
}

void SPH_main::place_points_beach(double* min, double* max)
{
	double x[2] = { min[0], min[1] };
	SPH_particle particle;
	double temp_pp;
	// place boundary particles
	// array that holds the range x-values for the leftmost points of the boundary
	double init_x_b[5] = { 0, (max[0] - 2 * h) * (7.0 / 8.0), min[0], max[0] - 2.0 * h, 0, };
	// array that holds the range of y-values for the leftmost points of the boundary
	double init_y_b[5] = { min[1], (max[0] - 2 * h) * (1.0 / 8.0) - 2 * h, min[0], (max[0] - 2 * h) * (1.0 / 8.0) - 2 * h, max[1] - 1.7 * h, };
	// array that holds the range x-values for the rightmost points of the boundary
	double set_x_b[5] = { (max[0] - 2 * h) * 3.0 / 4.0, max[0] - 2 * h, 0, max[0], max[0] - 2.0 * h , };
	// array that holds the range y-values for the rightmost points of the boundary
	double set_y_b[5] = { 0, (max[0] - 2 * h) * (1.0 / 8.0), max[1],  max[1], max[1], };
	for (int i = 0; i < 5; i++)
	{
		x[0] = init_x_b[i];
		x[1] = init_y_b[i];
		while (x[0] < set_x_b[i])
		{
			x[1] = init_y_b[i];
			while (x[1] < set_y_b[i])
			{
				for (int i = 0; i < 2; i++)
					particle.x[i] = x[i];
				particle.boundary = true;
				particle.calc_index();
				particle_list.push_back(particle);
				x[1] += dx;
			}
			x[0] += dx;
		}
	}
	// bottom-slope
	x[0] = (max[0] - 2 * h) * (3.0 / 4.0);
	temp_pp = min[1];
	while (x[0] <= ((max[0] - 2 * h) * (7.0 / 8.0)))
	{
		x[1] = temp_pp;
		while (x[1] <= (temp_pp + 2 * h))
		{
			for (int i = 0; i < 2; i++)
				particle.x[i] = x[i];
			particle.calc_index();
			particle.boundary = true;
			particle_list.push_back(particle);
			x[1] += dx;
		}
		temp_pp += dx;
		x[0] += dx;
	}
	// place fluid particle
	// array that holds the range x-values for the leftmost points of the lake
	double init_x[2] = { dx / 2, dx / 2 };
	// array that holds the range of y-values for the leftmost points of the lake
	double init_y[2] = { dx / 2, 2.0 + dx / 2 };
	// array that holds the range x-values for the rightmost points of the lake
	double set_x[2] = { (max[0] - 2.0 * h) * 3.0 / 4.0, 3 };
	// array that holds the range y-values for the rightmost points of the lake
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
					particle.x[i] = x[i];
				particle.boundary = false;
				particle.calc_index();
				particle_list.push_back(particle);
				x[1] += dx;
			}
			x[0] += dx;
		}
	}
	// slope fluid particle
	x[0] = (max[0] - 2 * h) * (3.0 / 4.0);// +dx / 2;
	temp_pp = dx / 2;
	while (x[0] <= ((max[0] - 2 * h) * (3.0 / 4.0) + set_y[0]))
	{
		x[1] = temp_pp;
		while (x[1] <= set_y[0])
		{
			for (int i = 0; i < 2; i++)
				particle.x[i] = x[i];
			particle.calc_index();
			particle.boundary = false;
			particle_list.push_back(particle);
			x[1] += dx;
		}
		temp_pp += dx;
		x[0] += dx;
	}
}



void SPH_main::allocate_to_grid(void)				//needs to be called each time that all the particles have their positions updated
{
	for (int i = 0; i < max_list[0]; i++)
		for (int j = 0; j < max_list[1]; j++)
			search_grid[i][j].clear();

	for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
	{
		particle_list[cnt].part_num = search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].size();
		search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
	}
}

void SPH_main::neighbour_iterate(SPH_particle* part)					//iterates over all particles within 2h of part - can be made more efficient using a stencil and realising that all interactions are symmetric
{
	SPH_particle* other_part;
	double dist;			//distance between particles
	double dn[2];			//vector from 1st to 2nd particle
	double tmp_CFL_dt;
	double e_ij[2], v_ij[2];

	double f_rep = 0;

	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0]) {
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1])
				{
					// calculate dist_to_edge x and y for particle

					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
					{
						other_part = search_grid[i][j][cnt];	// re-direct the address

						if (part != other_part)					//stops particle interacting with itself
						{
							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
								dn[n] = part->x[n] - other_part->x[n];

							dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

							if (dist < 2. * h && dist > 0.)					//only particle within 2h
							{
								//TODO: all the interactions between the particles
								for (int n = 0; n < 2; n++)
								{
									v_ij[n] = part->v[n] - other_part->v[n];
									e_ij[n] = (part->x[n] - other_part->x[n]) / dist;
								}

								// Calculate the Courant number to check timestep is still valid
								tmp_CFL_dt = h / sqrt((v_ij[0] * v_ij[0]) + (v_ij[1] * v_ij[1]));
								if (tmp_CFL_dt < this->dt_CFL) {
									this->dt_CFL = tmp_CFL_dt;
								}

								part->NS_acceleration(other_part, e_ij, v_ij, dist);
								part->NS_density(other_part, e_ij, v_ij, dist);

							}
						}
					}
				}
		}

}

void SPH_main::neighbour_iterate_stencil(SPH_particle* part)					//iterates over all particles within 2h of part - can be made more efficient using a stencil and realising that all interactions are symmetric
{
	SPH_particle* other_part;
	double dist;			//distance between particles
	double dn[2];			//vector from 1st to 2nd particle

	signed int row_nei[5] = { 0, 1, -1, 0, 1 };
	signed int col_nei[5] = { 0, 0,  1, 1, 1 };


	double D = this->g * this->h;
	double n1 = 4;
	double n2 = 2;
	double cut_off = 0.5 * this->dx;
	double tmp_CFL_dt;

	double e_ij[2], v_ij[2];

	//cout << "part->list_num[0]: " << part->list_num[0] << " part->list_num[1]: " << part->list_num[1] << endl;
	for (int k = 0; k < 5; k++)
	{
		int i = part->list_num[0] + row_nei[k];
		int j = part->list_num[1] + col_nei[k];
		//cout << "i: " << i << " j: " << j << endl;
		if (i >= 0 && i < max_list[0] && j >= 0 && j < max_list[1])
		{
			if (k != 0)
			{
				for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
				{
					other_part = search_grid[i][j][cnt];

					//Calculates the distance between potential neighbours
					for (int n = 0; n < 2; n++)
						dn[n] = part->x[n] - other_part->x[n];

					dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

					if (dist < 2. * h && dist > 0.)					//only particle within 2h
					{
						// All the interactions between the particles
						for (int n = 0; n < 2; n++)
						{
							v_ij[n] = part->v[n] - other_part->v[n];
							e_ij[n] = (part->x[n] - other_part->x[n]) / dist;
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

					//Calculates the distance between potential neighbours
					for (int n = 0; n < 2; n++)
						dn[n] = part->x[n] - other_part->x[n];

					dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

					if (dist < 2. * h && dist > 0.)					//only particle within 2h
					{
						// All the interactions between the particles
						for (int n = 0; n < 2; n++)
						{
							v_ij[n] = part->v[n] - other_part->v[n];
							e_ij[n] = (part->x[n] - other_part->x[n]) / dist;
						}

						// Calculate the Courant number to check timestep is still valid
						tmp_CFL_dt = h / sqrt((v_ij[0] * v_ij[0]) + (v_ij[1] * v_ij[1]));
						if (tmp_CFL_dt < this->dt_CFL) {
							this->dt_CFL = tmp_CFL_dt;
						}
						part->NS_acceleration_simultaneous(other_part, e_ij, v_ij, dist);
						part->NS_density_simultaneous(other_part, e_ij, v_ij, dist);
					}

				}
			}
		}
	}

}
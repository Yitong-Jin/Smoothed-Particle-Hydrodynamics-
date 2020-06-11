#include "SPH_2D.h"
#include "file_writer.h"
#include <math.h>
#include <string>
#include <omp.h>

#define M_PI 3.14159265358979323846

using namespace std;

SPH_main domain;

// Set domain type (whether replicates a beach or not) - comment as appropriate
//string domain_type = "beach";
string domain_type = "no beach";

int main(void)
{
	domain.set_values();										//Set simulation parameters

	domain.initialise_grid();									//initialise simulation grid

	// Change domain instantiation dependent on what domain the user prefers
	if (domain_type == "beach") { domain.place_points_beach(domain.min_x, domain.max_x); }			//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
	else { domain.place_points(domain.min_x, domain.max_x); }

	// Write out initial domain to file
	const string fname = "../tests/step_" + to_string(0) + ".vtp";
	write_file(fname.c_str(), &domain.particle_list);

	double time = 0;
	double max_time = 30;
	double print_every = 0.1;				// Print-out interval in seconds
	int file_writer = 1;					// files written incrementor
	double check_point = print_every;		// checkpoint marker - will write to file whenever we first exceed that checkpoint
											// and then will increment the checkpoint
	int it = 0;								// while loop count

	// Instantiate the particle variables 
	for (int i = 0; i < domain.particle_list.size(); i++)
	{
		for (int n = 0; n < 2; n++)
		{
			domain.particle_list[i].x[n][1] = domain.particle_list[i].x[n][0];
			domain.particle_list[i].v[n][1] = 0.0;
			domain.particle_list[i].rho[1] = domain.particle_list[i].rho[0];
		}
	}

	// Initiate time measurements
	double t_start = omp_get_wtime();
	while (time < max_time)
	{
		//needs to be called for each time step
		for (int step = 1; step >= 0; step--)
		{
			// Allocate points to their grids in space
			domain.allocate_to_grid();

			// Smooth densities every 10-20 timesteps
			if (it % 20 == 0)
				domain.smoothing();


#pragma omp parallel for schedule(guided)										// deconstruct the domain into parallel regions with omp parallel for, telling scheduler that we want to guide the schedule with larger values first	
			for (int i = 0; i < domain.particle_list.size(); i++)
			{
				domain.neighbour_iterate_stencil(&domain.particle_list[i]);		//finding all the neighbours of the 100th particle in the list - in reality the simulation loop will need to do the calculations for the neighbours of every particle															
			}

#pragma omp parallel for schedule(guided)										// deconstruct the domain into parallel regions with omp parallel for, telling scheduler that we want to guide the schedule with larger values first
			for (int i = 0; i < domain.particle_list.size(); i++)
			{
				domain.particle_list[i].update_pc(step, domain_type);			// do predictor-corrector update to complete the timesteps
			}
		}

		if (time > check_point)													// if the current time has exceeded the checkpoint ...
		{
			cout << time << " for iteration " << it << endl;
			const string fname = "../tests/step_" + to_string(file_writer) + ".vtp";
			write_file(fname.c_str(), &domain.particle_list);					// ... write out to file
			check_point += print_every;											// set new checkpoint
			file_writer++;														// increment the files-written counter
		}

		time += domain.dt;														// update the time using the adaptive dt
		it++;
	}
	double t_end = omp_get_wtime();												// finish time measurements
	cout << "Time taken = " << t_end - t_start << endl;

	return 0;
}

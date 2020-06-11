#include "../includes/file_writer.h"
#include "../includes/SPH_2D.h"
#include <math.h>
#include <string>
#include <omp.h>

#define M_PI 3.14159265358979323846

using namespace std;
///string keyword = "beach";
string keyword = "no beach";
SPH_main domain;

int main(void)
{
	domain.set_values();										//Set simulation parameters

	domain.initialise_grid();									//initialise simulation grid

	if (keyword == "beach")
		domain.place_points_beach(domain.min_x, domain.max_x);			//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
	else
		domain.place_points(domain.min_x, domain.max_x);
	const string fname = "./Forward_Euler_repo/tests/step" + to_string(0) + ".vtp";
	write_file(fname.c_str(), &domain.particle_list);

	double t_start = omp_get_wtime();
	double time = 0;
	double max_time = 20;
	double print_every = 0.05;
	int file_writer = 1;
	double check_point = print_every;
	int it = 0;

	while (time < max_time)
	{
		domain.allocate_to_grid();			//needs to be called for each time step

		if (it % 20 == 0) {
			domain.smoothing();
		}
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < domain.particle_list.size(); i++)
		{
			domain.neighbour_iterate_stencil(&domain.particle_list[i]);		//finding all the neighbours of the 100th particle in the list - in reality the simulation loop will need to do the calculations for the neighbours of every particle															
		}

#pragma omp parallel for schedule(guided)
		for (int i = 0; i < domain.particle_list.size(); i++)
		{
			domain.particle_list[i].update(keyword);
		}

		if (time > check_point)
		{
			cout << time << endl;
			const string fname = "./Forward_Euler_repo/tests/step_" + to_string(file_writer) + ".vtp";
			write_file(fname.c_str(), &domain.particle_list);
			check_point += print_every;
			file_writer++;
		}

		time += domain.dt;
		it++;
	}
	double t_end = omp_get_wtime();
	cout << "Time taken = " << t_end - t_start << endl;


	return 0;
}
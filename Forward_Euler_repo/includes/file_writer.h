#pragma once

#include <vector>
#include "../includes/SPH_2D.h"

int write_file(const char* filename,
	std::vector<SPH_particle>* particle_list);
#ifndef SOLVER_HPP
#define SOLVER_HPP
#include <vector>
#include "model.hpp"

class solver_edp
{
public:
	
	solver_edp(model pde_model);
	std::vector<double> solve_pde();
private:
	model s_pde_model;
	
};

#endif
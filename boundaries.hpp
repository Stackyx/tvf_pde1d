#ifndef BOUND_HPP
#define BOUND_HPP
#include <vector>
#include "model.hpp"
#include "mesh.hpp"

class bound
{
public:

	bound(payoff f, mesh grille, std::string method = "Dirichlet", std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}});
	
	void get_boundaries(double r, double T, double t, double& b_down, double& b_up);
	
private:
	std::vector<std::vector<double>> bound::getNeumann();
	std::vector<std::vector<double>> bound::getDirichelet();
	std::vector<std::vector<double>> bound::get_conditions(std::vector<std::vector<double>> conditions, std::string method);
	
	payoff b_f;
	mesh b_mesh;
	
	std::string b_method;
	std::vector<double> strikes;
	
};
	

#endif
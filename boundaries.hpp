#ifndef BOUND_HPP
#define BOUND_HPP
#include <vector>
#include "model.hpp"
#include "mesh.hpp"

class bound
{
public:

	bound(model pde_model, mesh grille, std::string method = "Dirichlet", std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}});
	
	std::vector<std::vector<double>> get_boundaries();
	
private:
	std::vector<std::vector<double>> bound::getNeumann();
	std::vector<std::vector<double>> bound::getDirichelet();
	std::vector<std::vector<double>> bound::get_conditions(std::vector<std::vector<double>> conditions, std::string method);
	
	model m_pde_model;
	mesh m_grille;
	std::string m_method;
	std::vector<std::vector<double>> m_cdt;
};
	

#endif
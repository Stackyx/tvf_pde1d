#ifndef BOUND_HPP
#define BOUND_HPP
#include <vector>
#include "model.hpp"
#include "mesh.hpp"

class bound
{
public:

	bound(const payoff& f, const mesh& grille, const std::vector<double>& conditions, std::string method = "Dirichlet");
	bound(const payoff& f, const mesh& grille, std::string method = "Dirichlet");
	
	void adapt_mat(std::vector<std::vector<double>>& mat, std::vector<std::vector<double>>& mat_inv, double theta, double r, const std::vector<double>& sigma) const;
	void get_boundaries(std::vector<double>& sol, double T, double dt, int i, double r) const;
		
private:
	
	const payoff& b_f;
	const mesh& b_mesh;
	
	std::string b_method;
	double b_conditions_up;
	double b_conditions_down;
};
	

#endif
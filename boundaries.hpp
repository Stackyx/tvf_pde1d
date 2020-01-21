#ifndef BOUND_HPP
#define BOUND_HPP
#include <vector>
#include "model.hpp"
#include "mesh.hpp"

class bound
{
public:

	bound(payoff f, mesh grille, std::string method, std::vector<std::vector<double>> conditions);
	bound(payoff f, mesh grille, std::string method);
	bound(payoff f, mesh grille, std::vector<std::vector<double>> conditions);
	bound(payoff f, mesh grille);
	
	void adapt_mat(std::vector<std::vector<double>>& mat_inv);
	void get_boundaries(double ri, double ri1, double sigma0, double sigma1, double T, double dt, double j, std::vector<double>& sol, const std::vector<double>& sol_back);
	
private:
	void get_boundaries_cdt(double ri, double ri1, double sigma0, double sigma1, double T, double dt, double j, std::vector<double>& sol, const std::vector<double>& sol_back);
	void get_boundaries_nocdt(double ri, double ri1,double sigma0, double sigma1, double T, double	dt, double j, std::vector<double>& sol, const std::vector<double>& sol_back);
	
	payoff b_f;
	mesh b_mesh;
	
	std::string b_method;
	std::vector<double> b_strikes;
	std::vector<std::vector<double>> b_conditions;
	
};
	

#endif
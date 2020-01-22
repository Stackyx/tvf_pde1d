#ifndef BOUND_HPP
#define BOUND_HPP
#include <vector>
#include "model.hpp"
#include "mesh.hpp"

class bound
{
public:

	bound(payoff f, mesh grille, std::vector<double> conditions, std::string method = "Dirichlet");
	bound(payoff f, mesh grille, std::string method = "Dirichlet");
	
	void adapt_mat(std::vector<std::vector<double>>& mat, std::vector<std::vector<double>>& mat_inv, double theta, double r, std::vector<double> sigma);
	//void get_boundaries(double ri, double ri1, double sigma0, double sigma1, double T, double dt, double j, std::vector<double>& sol, const std::vector<double>& sol_back);
	void get_boundaries(std::vector<double>& sol, double T, double dt, int i, double r);
	
	double get_cdt_up();
	double get_cdt_down();
	
private:
	//void get_boundaries_cdt(double ri, double ri1, double sigma0, double sigma1, double T, double dt, double j, std::vector<double>& sol, const std::vector<double>& sol_back);
	//void get_boundaries_nocdt(double ri, double ri1,double sigma0, double sigma1, double T, double	dt, double j, std::vector<double>& sol, const std::vector<double>& sol_back);
	
	payoff b_f;
	mesh b_mesh;
	
	std::string b_method;
	double b_conditions_up;
	double b_conditions_down;
};
	

#endif
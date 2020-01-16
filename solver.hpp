#ifndef SOLVER_HPP
#define SOLVER_HPP
#include <vector>
#include "model.hpp"

class solver_edp
{
public:
	
	solver_edp(model pde_model, std::string method = "Dirichlet", std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}});
	
	void solve_pde();
		
	std::vector<double> solution;
	std::vector<double> delta;
	std::vector<double> gamma;

private:
	model s_pde_model;
	
	void pde_matrix(std::vector<std::vector<double>>& mat, const std::vector<double>& sigma, const double& r, const double& theta, const double& dt, const double& dx, const int& nx, const int& i);
	void pde_matrix_to_inverse(std::vector<std::vector<double>>& mat, const std::vector<double>& sigma, const double& r, const double& theta, const double& dt, const double& dx, const int& nx, const int& i);
	void product_inverse(std::vector<double>& x, std::vector<std::vector<double>> trig_mat, std::vector<double> d);
	void trig_matmul(std::vector<double>& res, std::vector<std::vector<double>> trig_mat, std::vector<double> x);
	
	std::string s_method;
	
	std::vector<std::vector<double>> s_cdt;
	std::vector<std::vector<double>> getDirichelet();
	std::vector<std::vector<double>> getNeumann();
	std::vector<std::vector<double>> get_conditions(std::vector<std::vector<double>> conditions, std::string method);

};

#endif
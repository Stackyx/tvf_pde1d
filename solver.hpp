#ifndef SOLVER_HPP
#define SOLVER_HPP
#include <vector>
#include "model.hpp"

class solver_edp
{
public:
	
	solver_edp(model pde_model);
	std::vector<double> solve_pde();
	
	void product_inverse(std::vector<double>& x, std::vector<std::vector<double>> trig_mat, std::vector<double> d);
	void trig_matmul(std::vector<double>& res, std::vector<std::vector<double>> trig_mat, std::vector<double> x);
private:
	model s_pde_model;
	void pde_matrix(std::vector<std::vector<double>>& mat, const double& sigma, const double& r, const double& theta, const double& dt, const double& dx, const int& nx, const int& i);
	void pde_matrix_to_inverse(std::vector<std::vector<double>>& mat, const double& sigma, const double& r, const double& theta, const double& dt, const double& dx, const int& nx, const int& i);
};

#endif
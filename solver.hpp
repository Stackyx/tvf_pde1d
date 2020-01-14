#ifndef SOLVER_HPP
#define SOLVER_HPP
#include <vector>
#include "model.hpp"

class solver_edp
{
public:
	
	solver_edp(model pde_model);
	std::vector<double> solve_pde();
	
	std::vector<double> product_inverse(std::vector<std::vector<double>> trig_mat, std::vector<double> d);
	std::vector<double> trig_matmul(const std::vector<std::vector<double>>& trig_mat, const std::vector<double>& x);
private:
	model s_pde_model;
	
};


class greeks
{
public:
	greeks(solver_edp rst_solver);
	
	std::vector<double> getDelta();
	std::vector<double> getGamma();
	std::vector<double> getTheta();
	std::vector<double> getVega();
	
private:
	solver_edp m_solver;
}
#endif
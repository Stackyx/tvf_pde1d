#ifndef SOLVER_HPP
#define SOLVER_HPP
#include <vector>
#include "model.hpp"
#include "mesh.hpp"
#include "boundaries.hpp"

class solver_edp
{
public:
	
	solver_edp(model pde_model, mesh grille, bound boundary, payoff f, double theta);
	
	void solve_pde(const bool& vega_bool = 0);
	void export_csv(std::string f_name = "output.csv");
	void print_results();
	
	std::vector<double> solution;
	std::vector<double> delta;
	std::vector<double> gamma;
	std::vector<double> vega;

private:
	model s_pde_model;
	mesh s_mesh;
	bound s_bound;
	payoff s_f;
	
	double s_theta;
	
	void pde_matrix(std::vector<std::vector<double>>& mat, std::vector<std::vector<double>>& mat_inv, const std::vector<double>& sigma, const std::vector<double>& sigma_plus, double r, double r_plus, double theta, double dt, double dx, int nx, int i);
	void product_inverse(std::vector<double>& x, std::vector<std::vector<double>> trig_mat, std::vector<double> d);
	void trig_matmul(std::vector<double>& res, std::vector<std::vector<double>> trig_mat, std::vector<double> x);
	
	std::vector<std::vector<double>> s_cdt;
	

};

#endif

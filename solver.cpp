#include <vector>
#include "solver.hpp"
#include <iostream>

solver_edp::solver_edp(model pde_model)
	: s_pde_model(pde_model)
{
}

std::vector<double> solver_edp::solve_pde()
{	
	std::vector<double> a;
	return a;
}

std::vector<double> solver_edp::product_inverse(std::vector<std::vector<double>> trig_mat, std::vector<double> d)
{
	double w;
	std::vector<double> a(trig_mat[0]);
	std::vector<double> b(trig_mat[1]);
	std::vector<double> c(trig_mat[2]);
	
	std::vector<double> x(trig_mat[0].size());
	
	for(int i=1; i< trig_mat[0].size(); ++i)
	{
		w = a[i]/b[i-1];
		b[i] = b[i] - w*c[i-1];
		d[i] = d[i] - w*d[i-1];
	}
	
	x[trig_mat[0].size()-1] = d[trig_mat[0].size()-1]/b[trig_mat[0].size()-1];

	for(int i=trig_mat[0].size()-2; i >= 0; --i)
	{
		x[i] = (d[i]-c[i]*x[i+1])/b[i];
	}
	
	return x;
}

std::vector<double> solver_edp::trig_matmul(const std::vector<std::vector<double>>& trig_mat, const std::vector<double>& x)
{
	std::vector<double> inf(trig_mat[0]);
	std::vector<double> diag(trig_mat[1]);
	std::vector<double> sup(trig_mat[2]);
	
	std::vector<double> res(trig_mat[0]);
	
	for(int i=1;i<inf.size()-2;i++)
	{
		res[i] = diag[i]*x[i] + inf[i]*x[i-1] + sup[i]*x[i+1];
	}
	
	res[0] = diag[0]*x[0] + sup[0]*x[1];
	res[res.size()-1] = diag[res.size()-1]*x[res.size()-1] + inf[res.size()-1]*x[res.size()-1];
	
	return res;
}
#include <vector>
#include "solver.hpp"
#include <iostream>

solver_edp::solver_edp(model pde_model)
	: s_pde_model(pde_model)
{
}

std::vector<double> solver_edp::solve_pde()
{	
	std::vector<std::vector<double>> boundaries(s_pde_model.m_cdt);
	
	std::vector<double> res(s_pde_model.m_nx);
	res[0] = boundaries[0][s_pde_model.m_nt-1];
	res[s_pde_model.m_nx-1] = boundaries[1][s_pde_model.m_nt-1];
	
	std::vector<std::vector<double>> pde_mat(3, std::vector<double>(s_pde_model.m_nx));
	std::vector<std::vector<double>> pde_mat_inv(3, std::vector<double>(s_pde_model.m_nx));
	
	for(int i=1; i< res.size()-1; ++i)
	{
		res[i] = s_pde_model.m_f.getpayoff()(exp(s_pde_model.m_Smin+i*s_pde_model.m_dx));
	}

	for(int i = s_pde_model.m_nt - 1; i > 0; --i)
	{
		s_pde_model.pde_matrix_to_inverse(pde_mat_inv, i);
		s_pde_model.pde_matrix(pde_mat, i);
		
		res = product_inverse(pde_mat_inv, trig_matmul(pde_mat, res)); 
		
		res[0] = boundaries[0][i-1];
		res[s_pde_model.m_nx-1] = boundaries[1][i-1];
	}

	return res;
}

std::vector<double> solver_edp::product_inverse(std::vector<std::vector<double>> trig_mat, std::vector<double> d)
{
	double w;
	
	std::vector<double> x(trig_mat[0].size());
	
	for(int i=1; i< trig_mat[0].size(); ++i)
	{
		w = trig_mat[0][i]/trig_mat[1][i-1];
		trig_mat[1][i] = trig_mat[1][i] - w*trig_mat[2][i-1];
		d[i] = d[i] - w*d[i-1];
	}
	
	x[trig_mat[0].size()-1] = d[trig_mat[0].size()-1]/trig_mat[1][trig_mat[0].size()-1];
	
	for(int i=trig_mat[0].size()-2; i >= 0; --i)
	{
		x[i] = (d[i]-trig_mat[2][i]*x[i+1])/trig_mat[1][i];
	}
	
	return x;
}

std::vector<double> solver_edp::trig_matmul(std::vector<std::vector<double>> trig_mat, std::vector<double> x)
{
	
	std::vector<double> res(trig_mat[0].size());
	
	for(int i=1;i<trig_mat[0].size()-1;i++)
	{
		res[i] = trig_mat[1][i]*x[i] + trig_mat[0][i]*x[i-1] + trig_mat[2][i]*x[i+1];
	}
	
	res[0] = trig_mat[1][0]*x[0] + trig_mat[2][0]*x[1];
	res[res.size()-1] = trig_mat[1][res.size()-1]*x[res.size()-1] + trig_mat[0][res.size()-1]*x[res.size()-2];
	
	return res;
}





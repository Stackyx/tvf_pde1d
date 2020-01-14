#include <vector>
#include "solver.hpp"
#include <iostream>

solver_edp::solver_edp(model pde_model)
	: s_pde_model(pde_model)
{
}

std::vector<double> solver_edp::solve_pde()
{	
	double tau = s_pde_model.m_T;
	std::vector<std::vector<double>> boundaries(s_pde_model.m_cdt);
	
	std::vector<double> res(s_pde_model.m_nx);
	res[0] = boundaries[0][s_pde_model.m_nt-1];
	res[s_pde_model.m_nx-1] = boundaries[1][s_pde_model.m_nt-1];
	
	for(int i=1; i< res.size()-1; ++i)
	{
		res[i] = s_pde_model.m_f.getpayoff()(exp(s_pde_model.m_Smin+i*s_pde_model.m_dx));
	}

	for(int i = s_pde_model.m_nt - 1; i > 0; --i)
	{
		res = product_inverse(s_pde_model.pde_matrix_to_inverse(i), trig_matmul(s_pde_model.pde_matrix(i), res)); 
		
		res[0] = boundaries[0][i-1];
		res[s_pde_model.m_nx-1] = boundaries[1][i-1];
	}

	return res;
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
	
	std::vector<double> res(trig_mat[0].size());
	
	for(int i=1;i<inf.size()-1;i++)
	{
		res[i] = diag[i]*x[i] + inf[i]*x[i-1] + sup[i]*x[i+1];
	}
	
	res[0] = diag[0]*x[0] + sup[0]*x[1];
	res[res.size()-1] = diag[res.size()-1]*x[res.size()-1] + inf[res.size()-1]*x[res.size()-2];
	
	return res;
}

greeks::greeks(solver_edp rst_solver)
	:m_solver(rst_solver)
{
	
}

std::vector<double> getDelta()
{
	return 
}


std::vector<double> getGamma()
{
	
	return 
}
std::vector<double> getTheta()
{
	return 
}

std::vector<double> getVega()
{
	double bump = 0.01;
	std::vector<double> bump_vol = m_solver.m_sigma + bump;
	
	model bumpvol_model = model(m_solver.m_initS, bump_vol, m_solver.m_r, m_solver.m_T, m_solver.m_nt, m_solver.m_nx, m_solver.m_theta, m_solver.m_f, m_solver.m_cdt, m_solver.m_method);
	return 
}





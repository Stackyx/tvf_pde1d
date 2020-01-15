#include <vector>
#include "solver.hpp"
#include <iostream>

solver_edp::solver_edp(model pde_model)
	: s_pde_model(pde_model)
{
}

void solver_edp::solve_pde()
{	
	std::vector<std::vector<double>> boundaries(s_pde_model.m_cdt);
	
	std::vector<double> sol(s_pde_model.m_nx);
	std::vector<double> vect(s_pde_model.m_nx);
	std::vector<double> sigma, sigma_plus;
	
	double r;
	
	sol[0] = boundaries[0][s_pde_model.m_nt-1];
	sol[s_pde_model.m_nx-1] = boundaries[1][s_pde_model.m_nt-1];
	
	std::vector<std::vector<double>> pde_mat(3, std::vector<double>(s_pde_model.m_nx));
	std::vector<std::vector<double>> pde_mat_inv(3, std::vector<double>(s_pde_model.m_nx));
	
	for(int i=1; i< sol.size()-1; ++i)
	{
		sol[i] = s_pde_model.m_f.getpayoff()(exp(s_pde_model.m_Smin+i*s_pde_model.m_dx));
	}

	for(int i = s_pde_model.m_nt - 1; i > 0; --i)
	{
		
		sigma_plus = s_pde_model.get_vol_col(i);
		sigma = s_pde_model.get_vol_col(i-1);
		r = s_pde_model.get_r(i);
		
		pde_matrix_to_inverse(pde_mat_inv, sigma, r, s_pde_model.m_theta, s_pde_model.m_dt, s_pde_model.m_dx, s_pde_model.m_nx, i);
		pde_matrix(pde_mat, sigma_plus, r, s_pde_model.m_theta, s_pde_model.m_dt, s_pde_model.m_dx, s_pde_model.m_nx, i);
		
		trig_matmul(vect, pde_mat, sol);
		product_inverse(sol, pde_mat_inv, vect); 
		
		sol[0] = boundaries[0][i-1];
		sol[s_pde_model.m_nx-1] = boundaries[1][i-1];
	}
	
	delta.resize(sol.size()-2);
	gamma.resize(sol.size()-2);
	
	for(int i=1; i < sol.size()-1; ++i)
	{
		double dx = exp(s_pde_model.m_Smin+(i)*s_pde_model.m_dx) - exp(s_pde_model.m_Smin+(i-1)*s_pde_model.m_dx);
		delta[i-1] = (sol[i+1] - sol[i-1])/(2*dx);
		gamma[i-1] = (sol[i+1] - 2*sol[i] + sol[i-1])/pow(dx, 2);
	}
	
	solution = sol;
}

void solver_edp::pde_matrix(std::vector<std::vector<double>>& mat, const std::vector<double>& sigma, const double& r, const double& theta, const double& dt, const double& dx, const int& nx, const int& i)
{	
	mat[1][0] = 1;
	mat[1][nx-1] = 1;

	for(int j = 1; j<nx-1; ++j)
	{
		mat[1][j] = 1 - dt*(1-theta)*(pow(sigma[j]/dx,2) + r);
		mat[0][j] = dt*(1-theta)/(2*dx)*(pow(sigma[j],2)/dx + pow(sigma[j],2)/2. - r);
		mat[2][j] = dt*(1-theta)/(2*dx)*(pow(sigma[j],2)/dx - pow(sigma[j],2)/2. + r);
	}

}

void solver_edp::pde_matrix_to_inverse(std::vector<std::vector<double>>& mat, const std::vector<double>& sigma, const double& r, const double& theta, const double& dt, const double& dx, const int& nx, const int& i)
{
	mat[1][0] = 1;
	mat[1][nx-1] = 1;
	
	for(int j = 1; j<nx-1; ++j)
	{		
		mat[1][j] = 1+dt*theta*(pow(sigma[j]/dx,2) + r);
		mat[0][j] = dt*theta/(2*dx)*(-pow(sigma[j],2)/dx - pow(sigma[j],2)/2. + r);
		mat[2][j] = dt*theta/(2*dx)*(-pow(sigma[j],2)/dx + pow(sigma[j],2)/2. - r);
	}

}

void solver_edp::product_inverse(std::vector<double>& x, std::vector<std::vector<double>> trig_mat, std::vector<double> d)
{
	double w;
	
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
}

void solver_edp::trig_matmul(std::vector<double>& res, std::vector<std::vector<double>> trig_mat, std::vector<double> x)
{	
	for(int i=1;i<trig_mat[0].size()-1;i++)
	{
		res[i] = trig_mat[1][i]*x[i] + trig_mat[0][i]*x[i-1] + trig_mat[2][i]*x[i+1];
	}
	
	res[0] = trig_mat[1][0]*x[0] + trig_mat[2][0]*x[1];
	res[res.size()-1] = trig_mat[1][res.size()-1]*x[res.size()-1] + trig_mat[0][res.size()-1]*x[res.size()-2];
}





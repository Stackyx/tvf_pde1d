#include <vector>
#include "solver.hpp"
#include <iostream>
#include <algorithm>

solver_edp::solver_edp(model pde_model, mesh grille, bound boundary, payoff m, double theta)
	: s_pde_model(pde_model), s_bound(boundary), s_mesh(grille), s_f(m), s_theta(theta)
{
}

void solver_edp::solve_pde(const bool& vega_bool)
{	

	double dx(s_mesh.get_dx());
	double dt(s_mesh.get_dt());
	double s_min(s_mesh.get_Smin());
	double s_max(s_mesh.get_Smax());
	
	std::vector<double> sol(s_mesh.get_nx());
	std::vector<double> vect(s_mesh.get_nx());
	std::vector<double> sigma, sigma_plus;
	
	double r;
		
	s_bound.get_boundaries(s_pde_model.get_r(s_mesh.get_nt()-1), (s_mesh.get_nt()-1)*dt, (s_mesh.get_nt()-1)*dt, sol[0], sol[s_mesh.get_nx()-1]);
	
	std::vector<std::vector<double>> pde_mat(3, std::vector<double>(s_mesh.get_nx()));
	std::vector<std::vector<double>> pde_mat_inv(3, std::vector<double>(s_mesh.get_nx()));
		
	for(int i=1; i< sol.size()-1; ++i)
	{
		sol[i] = s_f.getpayoff()(exp(s_min + i*dx));
	}

	for(int i = s_mesh.get_nt() - 1; i > 0; --i)
	{
		
		sigma_plus = s_pde_model.get_vol_col(i);
		sigma = s_pde_model.get_vol_col(i-1);
		r = s_pde_model.get_r(i);
		
		pde_matrix_to_inverse(pde_mat_inv, sigma, r, s_theta, dt, dx, s_mesh.get_nx(), i);
		pde_matrix(pde_mat, sigma_plus, r, s_theta, dt, dx, s_mesh.get_nx(), i);
		
		trig_matmul(vect, pde_mat, sol);
		product_inverse(sol, pde_mat_inv, vect); 

		s_bound.get_boundaries(s_pde_model.get_r(i-1), (s_mesh.get_nt()-1)*dt, (i-1)*dt, sol[0], sol[s_mesh.get_nx()-1]);
		
		// Mettre neumann dans boundaries
		
		// else if (s_method == "Neumann")
		// {
			// sol[0] = (sol[0] - dt*(1./2.*sigma[0]*sigma[0]-r)*boundaries[0][i-1] + dt/2.*sigma[0]*sigma[0]*(2*sol[1]-2*boundaries[0][i-1]*dx)/(dx*dx))/(dt*sigma[0]*sigma[0]/(dx*dx)+r*dt+1);
			// sol[s_mesh.get_nx()-1] = (sol[s_mesh.get_nx()-1] - dt*(1./2.*sigma[s_mesh.get_nx()-1]*sigma[s_mesh.get_nx()-1]-r)*boundaries[1][i-1] + dt/2.*sigma[s_mesh.get_nx()-1]*sigma[s_mesh.get_nx()-1]*(2*sol[s_mesh.get_nx()-2]+2*boundaries[1][i-1]*dx)/(dx*dx))/(dt*sigma[s_mesh.get_nx()-1]*sigma[s_mesh.get_nx()-1]/(dx*dx)+r*dt+1);
		// }
		
	}
	
	delta.resize(sol.size()-2);
	gamma.resize(sol.size()-2);
	
	for(int i=1; i < sol.size()-2; ++i)
	{
		double dxi = exp(s_min+(i)*s_mesh.get_nx()) - exp(s_min + (i-1)*s_mesh.get_nx());
		double dxi1 = exp(s_min + (i+1)*s_mesh.get_nx()) - exp(s_min + i*s_mesh.get_nx());
		delta[i-1] = (sol[i+1] - sol[i-1])/(dxi+dxi1);
		gamma[i-1] = (sol[i+1] - 2*sol[i] + sol[i-1])/(dxi*dxi1);
	}
	
	solution = sol;
	
	// A corriger avec les nouvelles classes : 
	
	// if(vega_bool){
		// std::vector<double> sol2(s_mesh.get_nx());
		// std::vector<std::vector<double>> vol_bump = s_pde_model.m_sigma;
		// double h = 0.01;
		// auto lambda = [h](double d1) { return d1 + h; };

		// for(int c = 0; c<s_pde_model.m_sigma[0].size(); ++c)
		// {
			// std::transform(vol_bump[c].begin(), vol_bump[c].end(), vol_bump[c].begin(), lambda);
		// }
		// model model_bump_vol(s_pde_model.m_initS, vol_bump, s_pde_model.m_Smin, s_pde_model.m_Smax ,s_pde_model.m_r, s_pde_model.m_T, s_pde_model.m_nt , s_pde_model.m_nx, s_pde_model.m_theta, s_pde_model.m_f);
		// solver_edp sol2_edp(model_bump_vol, s_grille, s_method, s_cdt);
		// sol2_edp.solve_pde();
		// sol2 = sol2_edp.solution;

		// vega.resize(sol2.size());

		// for (int c = 0; c<sol2.size(); ++c)
		// {
			// vega[c] = (sol2[c] - solution[c])/h;
		// }

	// }
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

void solver_edp::pde_matrix_to_inverse(std::vector<std::vector<double>>& mat, const std::vector<double>& sigma, const double& r, const double& theta, const double& dt, const double& dx, const int& nx, const int& i){
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



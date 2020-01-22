#include <vector>
#include "solver.hpp"
#include "closed_form.hpp"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>

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
	
	double T = dt*(s_mesh.get_nt()-1);
	
	solution.resize(s_mesh.get_nx());
	
	std::vector<double> vect(s_mesh.get_nx());
	std::vector<double> sol_inter(s_mesh.get_nx());
	
	std::vector<double> sigma(s_mesh.get_nx()), sigma_plus(s_mesh.get_nx());
	
	double r(s_pde_model.get_r(s_mesh.get_nt()-1));
	double r_plus;
	
	s_pde_model.get_vol_col(sigma,s_mesh.get_nt()-1);
	
	s_bound.get_boundaries(solution, T, dt, s_mesh.get_nt()-1, r); // Compute terminal conditions
	
	std::vector<std::vector<double>> pde_mat(3, std::vector<double>(s_mesh.get_nx()));
	std::vector<std::vector<double>> pde_mat_inv(3, std::vector<double>(s_mesh.get_nx()));

	for(int i = s_mesh.get_nt() - 1; i > 0; --i)
	{
		
		s_pde_model.get_vol_col(sigma_plus, i);
		s_pde_model.get_vol_col(sigma, i-1);
		r_plus = s_pde_model.get_r(i);
		r = s_pde_model.get_r(i-1);
		
		pde_matrix(pde_mat, pde_mat_inv, sigma, sigma_plus, r, r_plus, s_theta, dt, dx, s_mesh.get_nx(), i);
		s_bound.adapt_mat(pde_mat,pde_mat_inv, s_theta, r, sigma); // Adapt the pde matrices to the boundary conditions
		
		trig_matmul(vect, pde_mat, solution);
		product_inverse(solution, pde_mat_inv, vect); 
		
		s_bound.get_boundaries(solution, T, dt, i-1, r);

	}
	
	delta.resize(solution.size());
	gamma.resize(solution.size());
	
	for(int i=1; i < solution.size()-1; ++i)
	{
		double dxi = exp(s_min+(i)*s_mesh.get_dx()) - exp(s_min + (i-1)*s_mesh.get_dx());
		double dxi1 = exp(s_min + (i+1)*s_mesh.get_dx()) - exp(s_min + i*s_mesh.get_dx());
		double dxi2 = exp(s_min+(i+1)*s_mesh.get_dx()) - exp(s_min+(i-1)*s_mesh.get_dx());
		delta[i] = (solution[i+1] - solution[i-1])/(dxi2);
		delta[0] = delta[1];
		//gamma[i] = (solution[i+1] - 2*solution[i] + solution[i-1])/(dxi*dxi);
		gamma[i] = (delta[i] - delta[i-1])/dxi;
		gamma[0] = gamma[1];
	}
	
	delta[solution.size()-1] = delta[solution.size()-2];
	gamma[solution.size()-1] = gamma[solution.size()-2];
	
	if(vega_bool){
		std::vector<double> sol2(s_mesh.get_nx());
		std::vector<std::vector<double>> vol_bump = s_pde_model.getSigma();
		double h = 0.01;
		auto lambda = [h](double d1) { return d1 + h; };

		for(int c = 0; c<vol_bump.size(); ++c)
		{
			std::transform(vol_bump[c].begin(), vol_bump[c].end(), vol_bump[c].begin(), lambda);
		}
		
		model model_bump_vol(vol_bump, s_pde_model.get_r());
		solver_edp sol2_edp(model_bump_vol, s_mesh, s_bound, s_f, s_theta);
		sol2_edp.solve_pde();
		sol2 = sol2_edp.solution;

		vega.resize(sol2.size());

		for (int c = 0; c<sol2.size(); ++c)
		{
			vega[c] = (sol2[c] - solution[c])/h;
		}

	}
}

void solver_edp::pde_matrix(std::vector<std::vector<double>>& mat, std::vector<std::vector<double>>& mat_inv, const std::vector<double>& sigma, const std::vector<double>& sigma_plus, double r, double r_plus, double theta, double dt, double dx, int nx, int i)
{	
	mat[1][0] = 1;
	mat[1][nx-1] = 1;

	mat_inv[1][0] = 1;
	mat_inv[1][nx-1] = 1;

	for(int j = 1; j<nx-1; ++j)
	{
		mat[1][j] = 1 - dt*(1-theta)*(pow(sigma_plus[j]/dx,2) + r_plus);
		mat[0][j] = dt*(1-theta)/(2*dx)*(pow(sigma_plus[j],2)/dx + pow(sigma_plus[j],2)/2. - r_plus);
		mat[2][j] = dt*(1-theta)/(2*dx)*(pow(sigma_plus[j],2)/dx - pow(sigma_plus[j],2)/2. + r_plus);
		
		mat_inv[1][j] = 1+dt*theta*(pow(sigma[j]/dx,2) + r);
		mat_inv[0][j] = dt*theta/(2*dx)*(-pow(sigma[j],2)/dx - pow(sigma[j],2)/2. + r);
		mat_inv[2][j] = dt*theta/(2*dx)*(-pow(sigma[j],2)/dx + pow(sigma[j],2)/2. - r);
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

void solver_edp::print_results()
{
	double sMin = s_mesh.get_Smin();
	double dx = s_mesh.get_dx();
	double prix;
	
	std::cout << std::fixed << std::setprecision( 5 );
	
	if (vega.empty())
	{
		for(int i=0; i<solution.size(); ++i)
		{
			prix = bs_price(exp(sMin+i*dx)/exp(-s_pde_model.get_r_avg()*s_mesh.get_mat()), s_mesh.get_S(), s_mesh.get_sigma(), s_mesh.get_mat(), 1)*exp(-s_pde_model.get_r_avg()*s_mesh.get_mat());
			std::cout << exp(sMin+i*dx) << ", sol = " << solution[i] << ", theory = " << prix << ", difference = " << prix - solution[i] << ", delta = " << delta[i] << ", gamma = " << gamma[i] << std::endl;
		}
	}
	else
	{
		for(int i=0; i<solution.size(); ++i)
		{
			prix = bs_price(exp(sMin+i*dx)/exp(-s_pde_model.get_r_avg()*s_mesh.get_mat()), s_mesh.get_S(), s_mesh.get_sigma(), s_mesh.get_mat(), 1)*exp(-s_pde_model.get_r_avg()*s_mesh.get_mat());
			std::cout << exp(sMin+i*dx) << ", sol = " << solution[i] << ", theory = " << prix << ", difference = " << prix - solution[i] << ", delta = " << delta[i] << ", gamma = " << gamma[i] << ", vega = " << vega[i] << std::endl;
		}
	}
		
}

void solver_edp::export_csv(std::string f_name)
{
	double sMin = s_mesh.get_Smin();
	double dx = s_mesh.get_dx();
	
	std::ofstream f(f_name);
	
	if (vega.empty())
	{
		f << "Spot,Price,Delta,Gamma" << "\n";
		
		for(int i = 0; i<solution.size(); ++i)
		{
			f <<  exp(sMin+i*dx) << "," << solution[i] << "," << delta[i] << "," << gamma[i] << "\n";
		}
	}
	else
	{
		f << "Spot,Price,Delta,Gamma,Vega" << "\n";
	
		for(int i = 0; i<solution.size(); ++i)
		{
			f <<  exp(sMin+i*dx) << "," << solution[i] << "," << delta[i] << "," << gamma[i] << "," << vega[i] << "\n";
		}
	}
	
	f.close();
}
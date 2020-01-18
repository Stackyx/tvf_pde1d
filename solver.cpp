#include <vector>
#include "solver.hpp"
#include <iostream>

solver_edp::solver_edp(model pde_model, mesh grille, std::string method, std::vector<std::vector<double>> conditions)
	: s_pde_model(pde_model), s_method(method)
{
	s_cdt = get_conditions(conditions, method);
}

void solver_edp::solve_pde(const bool& vega_bool)
{	
	std::vector<std::vector<double>> boundaries(s_cdt);
	
	std::vector<double> sol(s_pde_model.m_nx);
	std::vector<double> vect(s_pde_model.m_nx);
	std::vector<double> sigma, sigma_plus;
	double dt(s_pde_model.m_dt);
	double dx(s_pde_model.m_dx);
	double r;
		
	if (s_method == "Dirichlet")
	{
		sol[0] = boundaries[0][s_pde_model.m_nt-1];
		sol[s_pde_model.m_nx-1] = boundaries[1][s_pde_model.m_nt-1];
	}
	else if (s_method == "Neumann")
	{
		sol[0] = s_pde_model.getpayoff()(exp(s_pde_model.m_Smin));
		sol[s_pde_model.m_nx-1] = s_pde_model.getpayoff()(exp(s_pde_model.m_Smax));
	}
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
		
		if (s_method == "Dirichlet")
		{
			sol[0] = boundaries[0][i-1];
			sol[s_pde_model.m_nx-1] = boundaries[1][i-1];
		}
		else if (s_method == "Neumann")
		{
			sol[0] = (sol[0] - dt*(1./2.*sigma[0]*sigma[0]-r)*boundaries[0][i-1] + dt/2.*sigma[0]*sigma[0]*(2*sol[1]-2*boundaries[0][i-1]*dx)/(dx*dx))/(dt*sigma[0]*sigma[0]/(dx*dx)+r*dt+1);
			sol[s_pde_model.m_nx-1] = (sol[s_pde_model.m_nx-1] - dt*(1./2.*sigma[s_pde_model.m_nx-1]*sigma[s_pde_model.m_nx-1]-r)*boundaries[1][i-1] + dt/2.*sigma[s_pde_model.m_nx-1]*sigma[s_pde_model.m_nx-1]*(2*sol[s_pde_model.m_nx-2]+2*boundaries[1][i-1]*dx)/(dx*dx))/(dt*sigma[s_pde_model.m_nx-1]*sigma[s_pde_model.m_nx-1]/(dx*dx)+r*dt+1);
		}
		
	}
	
	delta.resize(sol.size()-2);
	gamma.resize(sol.size()-2);
	
	for(int i=1; i < sol.size()-2; ++i)
	{
		double dxi = exp(s_pde_model.m_Smin+(i)*s_pde_model.m_dx) - exp(s_pde_model.m_Smin+(i-1)*s_pde_model.m_dx);
		double dxi1 = exp(s_pde_model.m_Smin+(i+1)*s_pde_model.m_dx) - exp(s_pde_model.m_Smin+(i)*s_pde_model.m_dx);
		delta[i-1] = (sol[i+1] - sol[i-1])/(dxi+dxi1);
		gamma[i-1] = (sol[i+1] - 2*sol[i] + sol[i-1])/(dxi*dxi1);
	}
	
	solution = sol;
	
	if(vega_bool){
		std::vector<double> sol2(s_pde_model.m_nx);
		std::vector<std::vector<double>> vol_bump = s_pde_model.m_sigma;
		double h = 0.01;
		auto lambda = [](double d1) { return d1 + h; };

		for(int c = 0; c<s_pde_model.m_sigma[0].size(); ++c)
		{
			std::transform(vol_bump[c].begin(), vol_bump[c].end(), vol_bump[c].begin(), lambda);
		}
		model model_bump_vol = model_pde(s_pde_model.m_initS, vol_bump, s_pde_model.m_Smin, s_pde_model.m_Smax ,s_pde_model.m_r, s_pde_model.m_T, s_pde_model.m_nt , s_pde_model.m_nx, s_pde_model.m_theta, s_pde_model.m_f);
		sol2 = solver_edp(s_pde_model, s_method, s_cdt).solve_pde();

		vega.resize(sol2.size());

		for (int c = 0; c<solution2.size(); ++c)
		{
			vega[c] = (sol2 - solution)/h;
		}

	}
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

std::vector<std::vector<double>> solver_edp::getDirichelet()
{
	std::vector<double> cdt = s_pde_model.getStrike();
	std::vector<std::vector<double>> new_cdt;
	std::vector<double> uppercdt;
	std::vector<double> lowercdt;
	
	new_cdt.resize(s_pde_model.m_nt, std::vector<double>(cdt.size()));
	uppercdt.resize(s_pde_model.m_nt);
	lowercdt.resize(s_pde_model.m_nt);
	
	for (int j = 0; j<s_pde_model.m_nt; ++j)
	{
		for (int i = 0; i<cdt.size(); ++i)
		{
			new_cdt[j][i] = cdt[i] * exp(- s_pde_model.m_r[j] * s_pde_model.m_dt*j);
			//new_cdt[j][i] = cdt[i];


		}
		
		uppercdt[j] = payoff(s_pde_model.getName(), new_cdt[j]).getpayoff()(exp(s_pde_model.m_Smax));
		lowercdt[j] = payoff(s_pde_model.getName(), new_cdt[j]).getpayoff()(exp(s_pde_model.m_Smin));
		
	}
	
	return {lowercdt, uppercdt};
	
}


std::vector<std::vector<double>> solver_edp::getNeumann()
{

	double h = 0.000001;
	
	std::vector<double> cdt = s_pde_model.getStrike();
	std::vector<std::vector<double>> neu_cdt;
	std::vector<double> uppercdt;
	std::vector<double> lowercdt;
	
	neu_cdt.resize(s_pde_model.m_nt, std::vector<double>(cdt.size()));
	uppercdt.resize(s_pde_model.m_nt);
	lowercdt.resize(s_pde_model.m_nt);
	
	for (int j = 0; j<s_pde_model.m_nt; ++j)
	{
		for (int i = 0; i<cdt.size(); ++i)
		{
			neu_cdt[j][i] = cdt[i] * exp(- s_pde_model.m_r[j] * s_pde_model.m_dt*j);

		}
		
		uppercdt[j] = exp(s_pde_model.m_Smax)*((payoff(s_pde_model.getName(), neu_cdt[j])).getpayoff()(exp(s_pde_model.m_Smax) + h) - payoff(s_pde_model.getName(), neu_cdt[j]).getpayoff()(exp(s_pde_model.m_Smax)))/h;
		lowercdt[j] = exp(s_pde_model.m_Smin)*((payoff(s_pde_model.getName(), neu_cdt[j])).getpayoff()(exp(s_pde_model.m_Smin) + h) - payoff(s_pde_model.getName(), neu_cdt[j]).getpayoff()(exp(s_pde_model.m_Smin)))/h;

	}
	
	return {lowercdt, uppercdt};
	
}

std::vector<std::vector<double>> solver_edp::get_conditions(std::vector<std::vector<double>> conditions, std::string method)
{	
	
	std::vector<std::vector<double>> c = {{0, 0}, {0,0}};
	if (std::equal(conditions[0].begin(), conditions[0].end(), c[0].begin()) && std::equal(conditions[1].begin(), conditions[1].end(), c[1].begin()))
	{
		if (CaseSensitiveIsEqual(method, "Dirichlet"))
		{
			std::vector<std::vector<double>> drchlt = getDirichelet();
			conditions.resize(2, std::vector<double>(s_pde_model.m_nt));
			conditions = drchlt;
		}
		else if (CaseSensitiveIsEqual(method, "Neumann"))
		{
			std::vector<std::vector<double>> Neu = getNeumann();
			conditions.resize(2, std::vector<double>(s_pde_model.m_nt));
			conditions = Neu;
		}
		else
		{
			std::cout<< "Issue with the name of the method" << std::endl;
		}
	}
	
	return conditions;
}

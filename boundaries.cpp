#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "boundaries.hpp"

bound::bound(payoff f, mesh grille)
	: b_f(f), b_mesh(grille)
{
	 b_method = "Dirichlet";
	 b_conditions = {{0,0}, {0,0}};
}

bound::bound(payoff f, mesh grille, std::string method, std::vector<std::vector<double>> conditions)
	: b_f(f), b_mesh(grille), b_method(method), b_conditions(conditions)
{
	if ((!(CaseSensitiveIsEqual(b_method,"Dirichlet"))) && (!(CaseSensitiveIsEqual(b_method,"Neumann"))))
	{
		std::cout<< "please enter Dirichelet or Neumann only" << std::endl;
	}
}

bound::bound(payoff f, mesh grille, std::string method)
	: b_f(f), b_mesh(grille), b_method(method)
{
	b_conditions = {{0,0}, {0,0}};
	if ((!(CaseSensitiveIsEqual(b_method,"Dirichlet"))) && (!(CaseSensitiveIsEqual(b_method,"Neumann"))))
	{
		std::cout<< "please enter Dirichlet or Neumann only" << std::endl;
	}	
}

bound::bound(payoff f, mesh grille,std::vector<std::vector<double>> conditions)
	: b_f(f), b_mesh(grille), b_conditions(conditions)
{
	b_method = "Dirichlet";
}

void bound::adapt_mat(std::vector<std::vector<double>>& mat, std::vector<std::vector<double>>& mat_inv, double theta, double r, std::vector<double> sigma)
{
	if (CaseSensitiveIsEqual(b_method,"Neumann"))
	{
		
		double dx = b_mesh.get_dx();
		double dt = b_mesh.get_dt();
		
		mat[1][0] = -1/dx;
		mat[2][0] = 1/dx;
		mat[0][b_mesh.get_nx()-1] = -1/dx;
		mat[1][b_mesh.get_nx()-1] = 1/dx;
			
		mat_inv[1][1] = 1+dt*theta*(pow(sigma[1]/dx,2)*0.5 + r - (1./2.*sigma[1]*sigma[1] - r)/(2*dx));
		mat_inv[0][1] = -dt*theta/(2*dx)*(-pow(sigma[1],2)/dx - pow(sigma[1],2)/2. + r)*dx;
		
		mat_inv[1][b_mesh.get_nx()-2] = 1+dt*theta*(pow(sigma[b_mesh.get_nx()-2]/dx,2)*0.5 + r + (1./2.*sigma[b_mesh.get_nx()-2]*sigma[b_mesh.get_nx()-2] - r)/(2*dx));
		mat_inv[2][b_mesh.get_nx()-2] = dx * dt*theta/(2*dx)*(-pow(sigma[b_mesh.get_nx()-2],2)/dx + pow(sigma[b_mesh.get_nx()-2],2)/2. - r);
	}
	
	
}


void bound::get_boundaries(double ri, double ri1, double sigma0, double sigma1, double T, double dt, double j, std::vector<double>& sol, const std::vector<double>& sol_back)
{
	
	std::vector<std::vector<double>> c = {{0, 0}, {0,0}};
	if (std::equal(b_conditions[0].begin(), b_conditions[0].end(), c[0].begin()) && std::equal(b_conditions[1].begin(), b_conditions[1].end(), c[1].begin()))
	{	
		get_boundaries_nocdt(ri, ri1,sigma0, sigma1, T, dt, j, sol, sol_back);
	}
	else
	{
		get_boundaries_cdt(ri, ri1,sigma0, sigma1, T, dt, j, sol, sol_back);
	}
	
}

void bound::get_boundaries_nocdt(double ri, double ri1, double sigma0, double sigma1, double T, double dt, double j, std::vector<double>& sol, const std::vector<double>& sol_back)
{
	if (CaseSensitiveIsEqual(b_method,"Dirichlet"))
	{
		b_strikes = b_f.getparameters();
		
		for (int i = 0; i < b_strikes.size(); ++i)
		{
			b_strikes[i] = b_strikes[i] * std::exp(-ri1 * (T-dt*j));
		}
		
		sol[b_mesh.get_nx()-1] = payoff(b_f.getname(), b_strikes).getpayoff()(std::exp(b_mesh.get_Smax()));
		sol[0] = payoff(b_f.getname(), b_strikes).getpayoff()(std::exp(b_mesh.get_Smin()));
	}
	
	if (CaseSensitiveIsEqual(b_method,"Neumann"))
	{
		double h = 0.000001;
		b_strikes = b_f.getparameters();
		
		for (int i = 0; i < b_strikes.size(); ++i)
		{
			b_strikes[i] = b_strikes[i] * std::exp(-ri1 * (T-dt*j));
		}
		
		double temp_b_up = std::exp(b_mesh.get_Smax())*(payoff(b_f.getname(), b_strikes).getpayoff()(std::exp(b_mesh.get_Smax()) + h) - payoff(b_f.getname(), b_strikes).getpayoff()(std::exp(b_mesh.get_Smax())))/h;
		double temp_b_down = std::exp(b_mesh.get_Smin())*(payoff(b_f.getname(), b_strikes).getpayoff()(std::exp(b_mesh.get_Smin()) + h) - payoff(b_f.getname(), b_strikes).getpayoff()(std::exp(b_mesh.get_Smin())))/h;
		
		sol[0] = (sol_back[0] - dt*(1./2.*sigma0*sigma0-ri)*temp_b_down + dt/2.*sigma0*sigma0*(2*sol[1]-2*temp_b_down*b_mesh.get_dx())/(b_mesh.get_dx()*b_mesh.get_dx()))/(dt*sigma0*sigma0/(b_mesh.get_dx()*b_mesh.get_dx())+ri*dt+1);
		sol[b_mesh.get_nx()-1] = (sol_back[b_mesh.get_nx()-1] - dt*(1./2.*sigma1*sigma1-ri)*temp_b_up + dt/2.*sigma1*sigma1*(2*sol[b_mesh.get_nx()-2]+2*temp_b_up*b_mesh.get_dx())/(b_mesh.get_dx()*b_mesh.get_dx()))/(dt*sigma1*sigma1/(b_mesh.get_dx()*b_mesh.get_dx())+ri*dt+1);
	}
	
}

void bound::get_boundaries_cdt(double ri, double ri1, double sigma0, double sigma1, double T, double dt, double j, std::vector<double>& sol, const std::vector<double>& sol_back)
{
	
	if (CaseSensitiveIsEqual(b_method,"Dirichlet"))
	{
		sol[0] = b_conditions[j][0];
		sol[b_mesh.get_nx()-1] = b_conditions[j][1];
	
	}
	
	if (CaseSensitiveIsEqual(b_method,"Neumann"))
	{
			
		double temp_b_down = b_conditions[j][0];
		double temp_b_up = b_conditions[j][1];

		sol[0] = (sol_back[0] - dt*(1./2.*sigma0*sigma0-ri)*temp_b_down + dt/2.*sigma0*sigma0*(2*sol[1]-2*temp_b_down*b_mesh.get_dx())/(b_mesh.get_dx()*b_mesh.get_dx()))/(dt*sigma0*sigma0/(b_mesh.get_dx()*b_mesh.get_dx())+ri*dt+1);
		sol[b_mesh.get_nx()-1] = (sol_back[b_mesh.get_nx()-1] - dt*(1./2.*sigma1*sigma1-ri)*temp_b_up + dt/2.*sigma1*sigma1*(2*sol[b_mesh.get_nx()-2]+2*temp_b_up*b_mesh.get_dx())/(b_mesh.get_dx()*b_mesh.get_dx()))/(dt*sigma1*sigma1/(b_mesh.get_dx()*b_mesh.get_dx())+ri*dt+1);

	}
	
}

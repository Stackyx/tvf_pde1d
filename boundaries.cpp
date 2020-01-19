#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "boundaries.hpp"

bound::bound(payoff f, mesh grille)
	: b_f(f), b_mesh(grille)
{
	 b_method = "Dirichelet";
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
		std::cout<< "please enter Dirichelet or Neumann only" << std::endl;
	}	
}

bound::bound(payoff f, mesh grille,std::vector<std::vector<double>> conditions)
	: b_f(f), b_mesh(grille), b_conditions(conditions)
{
	b_method = "Dirichelet";
}



void bound::get_boundaries(double ri, double ri1, double sigma0, double sigma1, double T, double dt, double j, std::vector<double>& sol)
{
	
	std::vector<std::vector<double>> c = {{0, 0}, {0,0}};
	if (std::equal(b_conditions[0].begin(), b_conditions[0].end(), c[0].begin()) && std::equal(b_conditions[1].begin(), b_conditions[1].end(), c[1].begin()))
	{
		get_boundaries_nocdt(ri, ri1,sigma0, sigma1, T, dt, j, sol);
	}
	else
	{
		get_boundaries_cdt(j, sol);
	}
	
}

void bound::get_boundaries_nocdt(double ri, double ri1, double sigma0, double sigma1, double T, double dt, double j, std::vector<double>& sol)
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
		
		double temp_b_up = std::exp(b_mesh.get_Smin())*(payoff(b_f.getname(), b_strikes).getpayoff()(std::exp(b_mesh.get_Smax()) + h) - payoff(b_f.getname(), b_strikes).getpayoff()(std::exp(b_mesh.get_Smax())))/h;
		double temp_b_down = std::exp(b_mesh.get_Smin())*(payoff(b_f.getname(), b_strikes).getpayoff()(std::exp(b_mesh.get_Smin()) + h) - payoff(b_f.getname(), b_strikes).getpayoff()(std::exp(b_mesh.get_Smin())))/h;
		
		sol[0] = (sol[0] - dt*(1./2.*sigma0*sigma0-ri)*temp_b_down + dt/2.*sigma0*sigma0*(2*sol[1]-2*temp_b_down*b_mesh.get_nx())/(b_mesh.get_nx()*b_mesh.get_nx()))/(dt*sigma0*sigma0/(b_mesh.get_nx()*b_mesh.get_nx())+ri*dt+1);
		sol[b_mesh.get_nx()-1] = (sol[b_mesh.get_nx()-1] - dt*(1./2.*sigma1*sigma1-ri)*temp_b_up + dt/2.*sigma1*sigma1*(2*sol[b_mesh.get_nx()-2]+2*temp_b_up*b_mesh.get_nx())/(b_mesh.get_nx()*b_mesh.get_nx()))/(dt*sigma1*sigma1/(b_mesh.get_nx()*b_mesh.get_nx())+ri*dt+1);


	//comment avoir sigma

	}
	
}

void bound::get_boundaries_cdt(double j, std::vector<double>& sol)
{
	sol[0] = b_conditions[j][0];
	sol[b_mesh.get_nx()-1] = b_conditions[j][1];
	
}


// std::vector<std::vector<double>> bound::getDirichelet()
// {
	// std::vector<double> cdt = m_pde_model.getStrike();
	// std::vector<std::vector<double>> new_cdt;
	// std::vector<double> uppercdt;
	// std::vector<double> lowercdt;
	
	// new_cdt.resize(m_grille.m_nt, std::vector<double>(cdt.size()));
	// uppercdt.resize(m_grille.m_nt);
	// lowercdt.resize(m_grille.m_nt);
	
	// for (int j = 0; j<m_grille.m_nt; ++j)
	// {
		// for (int i = 0; i<cdt.size(); ++i)
		// {
			// new_cdt[j][i] = cdt[i] * std::exp(- m_pde_model.m_r[j] * m_grille.m_dt*j);

		// }
		
		// uppercdt[j] = payoff(m_pde_model.getName(), new_cdt[j]).getpayoff()(std::exp(m_grille.m_Smax));
		// lowercdt[j] = payoff(m_pde_model.getName(), new_cdt[j]).getpayoff()(std::exp(m_grille.m_Smin));
		
	// }
	
	// return {lowercdt, uppercdt};
	
// }


// std::vector<std::vector<double>> bound::getNeumann()
// {

	// double h = 0.000001;
	
	// std::vector<double> cdt = m_pde_model.getStrike();
	// std::vector<std::vector<double>> neu_cdt;
	// std::vector<double> uppercdt;
	// std::vector<double> lowercdt;
	
	// neu_cdt.resize(m_grille.m_nt, std::vector<double>(cdt.size()));
	// uppercdt.resize(m_grille.m_nt);
	// lowercdt.resize(m_grille.m_nt);
	
	// for (int j = 0; j<m_grille.m_nt; ++j)
	// {
		// for (int i = 0; i<cdt.size(); ++i)
		// {
			// neu_cdt[j][i] = cdt[i] * std::exp(- m_pde_model.m_r[j] * m_grille.m_dt*j);

		// }
		
		// uppercdt[j] = std::exp(m_grille.m_Smax)*((payoff(m_pde_model.getName(), neu_cdt[j])).getpayoff()(std::exp(m_grille.m_Smax) + h) - payoff(m_pde_model.getName(), neu_cdt[j]).getpayoff()(std::exp(m_grille.m_Smax)))/h;
		// lowercdt[j] = std::exp(m_grille.m_Smin)*((payoff(m_pde_model.getName(), neu_cdt[j])).getpayoff()(std::exp(m_grille.m_Smin) + h) - payoff(m_pde_model.getName(), neu_cdt[j]).getpayoff()(std::exp(m_grille.m_Smin)))/h;

	// }
	
	// return {lowercdt, uppercdt};
	
// }


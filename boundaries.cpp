#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "boundaries.hpp"

bound::bound(payoff f, mesh grille, std::string method, std::vector<std::vector<double>> conditions)
	: b_f(f), b_mesh(grille), b_method(method)
{
	//m_cdt = get_conditions(conditions, m_method);
}

void bound::get_boundaries(double r, double t, double dt, double& b_down, double& b_up)
{
	if (b_method == "Dirichlet")
	{
		std::vector<double> strike = b_f.getparameters();
	
		for (int i = 0; i < strike.size(); ++i)
		{
			strike[i] = strike[i] * std::exp(-r * t * dt);
		}
		
		b_up = payoff(b_f.getname(), strike).getpayoff()(std::exp(b_mesh.get_Smax()));
		b_down = payoff(b_f.getname(), strike).getpayoff()(std::exp(b_mesh.get_Smin()));
	}
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

// std::vector<std::vector<double>> bound::get_conditions(std::vector<std::vector<double>> conditions, std::string method)
// {	
	
	// std::vector<std::vector<double>> c = {{0, 0}, {0,0}};
	// if (std::equal(conditions[0].begin(), conditions[0].end(), c[0].begin()) && std::equal(conditions[1].begin(), conditions[1].end(), c[1].begin()))
	// {
		// if (CaseSensitiveIsEqual(method, "Dirichlet"))
		// {
			// std::vector<std::vector<double>> drchlt = getDirichelet();
			// conditions.resize(2, std::vector<double>(m_grille.m_nt));
			// conditions = drchlt;
		// }
		// else if (CaseSensitiveIsEqual(method, "Neumann"))
		// {
			// std::vector<std::vector<double>> Neu = getNeumann();
			// conditions.resize(2, std::vector<double>(m_grille.m_nt));
			// conditions = Neu;
		// }
		// else
		// {
			// std::cout<< "Issue with the name of the method" << std::endl;
		// }
	// }
	
	// return conditions;
// }
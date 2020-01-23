#include <iostream>
#include "closed_form.hpp"
#include "payoff.hpp"
#include "model.hpp"
#include "boundaries.hpp"
#include "solver.hpp"
#include "mesh.hpp"
#include <algorithm>

int main(int argc, char* argv[])
{	
	// Exemple d'utilisation de la class payoff

	//dauphine::payoff pp = dauphine::payoff("s", { 95 , 105 }, [](double d2) {return d2 * 2; }); //if the client wants to input his own function
	//std::function<double(double)> f = [&](double d1) { return std::max(d1 - 100, 0.0); };
	//payoff pp = payoff("Cal", {100},  f);
	payoff pp = payoff("Call", {100});
	std::cout << pp.getparameters()[0] << std::endl; // get the parameters
	
	double mat = 1;
	//double vol = 0.2;
	double S = 100;
	double theta = 1./2;
	//double r = 0.02;
	
	int nx = 1000;
	int nt = 252;
	
	std::vector<double> r(nt);
	std::vector<double> sigma(nt);
	
	for(int i=0; i<nt; ++i)
	{	
		//sigma[i] = std::max(i*.2/100.+0.2, 0.4);
		sigma[i] =0.2;
		//r[i]=std::max(i*.01/100., 0.04);
		r[i] = 0.02;
	}
	
	mesh grille(S, mat, nx, nt, 0.2);
	model model_pde(sigma, r, nt, nx);
	
	// std::vector<std::vector<double>> cdt(nt, std::vector<double>(2));
	
	// for (int i=0; i<nt; ++i)
	// {
		// cdt[i][0] = 0;
		// cdt[i][1] = 1.;
	// }
	
	bound boundaries(pp, grille, "Neumann");

	//model model_pde_r(100., 0.2, r, 1, 10, 10, 1./2, pp);

//	model model_pde_sigmar(100., sigma, r, 1, 10, 10, 1./2, pp);

	//model model_pde_sigma(100., sigma, 0.05, 1, 10, 10, 1./2, pp);
	
	solver_edp solver_model(model_pde, grille, boundaries, pp, theta);
	solver_model.solve_pde(1);
	//solver_model.export_csv();
	solver_model.print_results();
	
	return 0;

}

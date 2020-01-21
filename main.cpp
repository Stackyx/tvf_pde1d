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
	payoff pp = payoff("Call", {100}); //if we want to use the function already input in the class
	std::cout << pp.getparameters()[0] << std::endl; // get the parameters
	
	double mat = 1;
	double vol = 0.2;
	double S = 100;
	double theta = 1./2;
	double r = 0.00;
	
	int nx = 1000;
	int nt = 365;
	
	mesh grille(S, mat, nx, nt, vol);
	model model_pde(vol, r, nt, nx);
	
	// std::vector<std::vector<double>> cdt(nt, std::vector<double>(2));
	
	// for (int i=0; i<nt; ++i)
	// {
		// cdt[i][0] = 0;
		// cdt[i][1] = 1.;
	// }
	
	bound boundaries(pp, grille, "Neumann");
	
	// std::vector<double> r(10);
	// std::vector<double> sigma(10);
	
	// for(int i=0; i<10; ++i)
	// {	
		// sigma[i] = i*2./100.+0.2;
		// r[i]=i*3./100.;
	// }
	
	//model model_pde_r(100., 0.2, r, 1, 10, 10, 1./2, pp);

//	model model_pde_sigmar(100., sigma, r, 1, 10, 10, 1./2, pp);

	//model model_pde_sigma(100., sigma, 0.05, 1, 10, 10, 1./2, pp);
	
	solver_edp solver_model(model_pde, grille, boundaries, pp, theta);
	solver_model.solve_pde(1);
	
	// std::cout<< model_pde.getSigma().size() << std::endl;//taille colone
	// std::cout<< model_pde.getSigma()[0].size() << std::endl;//ligne taille
	// std::cout<< model_pde.get_vol_col(0).size() << std::endl;
	
	double sMin = grille.get_Smin();
	double dx = grille.get_dx();
	
	for(int i=0; i<solver_model.solution.size(); ++i)
	{
		double prix = bs_price(exp(sMin+i*dx)/exp(-r*mat), S, vol, mat, 1)*exp(-r*mat);
		std::cout << exp(sMin+i*dx) << ", sol = " << solver_model.solution[i] << ", theory = " << prix << ", difference = " << prix - solver_model.solution[i]<< ", vega = " << solver_model.vega[i]<<std::endl;
	}
	

	return 0;
}

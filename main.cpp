#include <iostream>
#include "closed_form.hpp"
#include "model.hpp"
#include "solver.hpp"
#include <algorithm>


int main(int argc, char* argv[])
{	
	// Exemple d'utilisation de la class payoff

	//dauphine::payoff pp = dauphine::payoff("s", { 95 , 105 }, [](double d2) {return d2 * 2; }); //if the client wants to input his own function
	payoff pp = payoff("call", {100}); //if we want to use the function already input in the class
	std::cout << pp.getparameters()[0] << std::endl; // get the parameters
	std::cout << pp.getname() << std::endl; // get the name
	std::cout << pp.getpayoff()(110) << std::endl; //get the payoff function and evaluate it at 90

	
	model model_pde(100., 0.2, 0.0, 1, 252, 200, 1./2., pp);
	
	std::vector<double> r(10);
	std::vector<double> sigma(10);
	
	for(int i=0; i<10; ++i)
	{	
		sigma[i] = i*2./100.+0.2;
		r[i]=i*3./100.;
	}
	
	//model model_pde_r(100., 0.2, r, 1, 10, 10, 1./2, pp);

//	model model_pde_sigmar(100., sigma, r, 1, 10, 10, 1./2, pp);

	//model model_pde_sigma(100., sigma, 0.05, 1, 10, 10, 1./2, pp);
	
	solver_edp solver_model(model_pde, "Neumann");
	solver_model.solve_pde();
	// std::cout<< model_pde.getSigma().size() << std::endl;//taille colone
	// std::cout<< model_pde.getSigma()[0].size() << std::endl;//ligne taille
	// std::cout<< model_pde.get_vol_col(0).size() << std::endl;
	
	double dx = model_pde.get_dx();
	double sMin = model_pde.getSmin();
	
	// std::cout << model_pde.get_dx() << std::endl;
	
	for(int i=0; i<solver_model.delta.size(); ++i)
	{
		std::cout << exp(sMin+i*dx) << ", " << solver_model.solution[i] << std::endl;
	}

	return 0;
}

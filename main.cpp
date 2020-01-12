#include <iostream>
#include "closed_form.hpp"
#include "model.hpp"
#include "solver.hpp"
#include <algorithm>


int main(int argc, char* argv[])
{	
	// Exemple d'utilisation de la class payoff

	//dauphine::payoff pp = dauphine::payoff("s", { 95 , 105 }, [](double d2) {return d2 * 2; }); //if the client wants to input his own function
	payoff pp = payoff("strangle", { 95 , 105 }); //if we want to use the function already input in the class
	std::cout << pp.getparameters()[0] << std::endl; // get the parameters
	std::cout << pp.getname() << std::endl; // get the name
	std::cout << pp.getpayoff()(90) << std::endl; //get the payoff function and evaluate it at 90
	std::cout << pp.getpayoff()(95) << std::endl;
	std::cout << pp.getpayoff()(100) << std::endl;
	std::cout << pp.getpayoff()(105) << std::endl;
	std::cout << pp.getpayoff()(115) << std::endl;
	
	model model_pde(10, 0.2, 0.05, 1./10., pp);
	solver_edp solver_model(model_pde);
	
	std::vector<double> sol(solver_model.solve_pde());
	
	return 0;
}

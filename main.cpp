#include <iostream>
#include "closed_form.hpp"
#include "payoff.hpp"
#include "model.hpp"
#include "boundaries.hpp"
#include "solver.hpp"
#include "mesh.hpp"
#include <algorithm>
#include <fstream>

int main(int argc, char* argv[])
{	
	// Exemple d'utilisation de la class payoff

	//dauphine::payoff pp = dauphine::payoff("s", { 95 , 105 }, [](double d2) {return d2 * 2; }); //if the client wants to input his own function
	//std::function<double(double)> f = [&](double d1) { return std::max(d1 - 100, 0.0); };
	//payoff pp = payoff("Cal", {100},  f);
	payoff pp = payoff("Call", {100});
	std::cout << pp.getparameters()[0] << std::endl; // get the parameters
	
	double mat = 10./252.;
	double vol = 0.2;
	double S = 100;
	double theta = 1./2;
	double r = 0.02;
	
	int nx = 2000;
	int nt = 10;
	
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
	
	double sMin = grille.get_Smin();
	double dx = grille.get_dx();
	
	// for(int i=0; i<solver_model.solution.size(); ++i)
	// {
		// double prix = bs_price(exp(sMin+i*dx)/exp(-r*mat), S, vol, mat, 1)*exp(-r*mat);
		// std::cout << exp(sMin+i*dx) << ", sol = " << solver_model.solution[i] << ", theory = " << prix << ", difference = " << prix - solver_model.solution[i]<< ", delta = " <<""<<std::endl;
	// }
	
	std::ofstream f("output.csv");
	// for(std::vector<double>::const_iterator i = solver_model.solution.begin(); i != solver_model.solution.end(); ++i) 
	// {
		// f << *i << '\n';
	// }
	
	f << "Spot,Price,Delta,Gamma,Vega" << "\n";
	
	for(int i = 0; i<solver_model.solution.size(); ++i)
	{
		f <<  exp(sMin+i*dx) << "," << solver_model.solution[i] << "," << solver_model.delta[i] << "," << solver_model.gamma[i] << "," << solver_model.vega[i] << "\n";
	}
	
	f.close();
	
	return 0;

}

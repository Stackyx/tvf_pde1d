#include <iostream>
#include "closed_form.hpp"
#include "payoff.hpp"
#include "model.hpp"
#include "boundaries.hpp"
#include "solver.hpp"
#include "mesh.hpp"
#include <algorithm>
#include <chrono>

int main(int argc, char* argv[])
{	

	//std::function<double(double)> f = [&](double d1) { return std::max(d1 - 100, 0.0); };

	payoff pp = payoff("Call", {100});
	
	double mat = 5;
	//double vol = 0.2;
	double S = 100;
	double theta = 1./2;
	//double r = 0.02;
	
	int nx = 1000;
	int nt = 252*5;
	
	std::vector<double> r(nt);
	std::vector<double> sigma(nt);
	
	for(int i=0; i<nt; ++i)
	{	
		//sigma[i] = std::max(i*.2/100.+0.2, 0.4);
		sigma[i] =0.2;
		//r[i]=std::max(i*.01/100., 0.04);
		r[i] = 0.02;
	}
	
	auto start = std::chrono::steady_clock::now();
	
	mesh grille(S, mat, nx, nt, 0.2);
	model model_pde(sigma, r, nt, nx);
	
	// std::vector<std::vector<double>> cdt(nt, std::vector<double>(2));
	
	// for (int i=0; i<nt; ++i)
	// {
		// cdt[i][0] = 0;
		// cdt[i][1] = 1.;
	// }
	
	bound boundaries(pp, grille, "Neumann");
	
	solver_edp solver_model(model_pde, grille, boundaries, pp, theta);
		
	solver_model.solve_pde(1);
	
	auto end = std::chrono::steady_clock::now();
	
	std::cout << "Time taken solving :" << std::chrono::duration <double, std::milli> (end - start).count() << " ms" << std::endl;
	
	//solver_model.export_csv();
	//solver_model.print_results();
	
	return 0;

}

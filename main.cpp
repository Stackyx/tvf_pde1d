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
	
	std::vector<double> cond{90, 110}; 
	
	model model_pde(100., 0.2, 0.05, 1, 1./10., pp, cond);
	
	std::cout << "ok" << std::endl;
	
	std::vector<double> r(10);
	std::vector<std::vector<double>> sigma(10, std::vector<double>(10));
	
	std::cout << "ok2" << std::endl;
	
	for(int i=0; i<10; ++i)
	{	
		for (int j=0; j<10; ++j)
		{
			sigma[i][j] = i*1./100. + 0.2;
		}
		r[i]=i*3./100.;
	}
	
	model model_pde_r(100., 0.2, r, 1, 1./10., pp, cond);
	std::cout << "ok3" << std::endl;
	model model_pde_sigmar(100., sigma, r, 1, 1./10., pp, cond);
	std::cout << "ok4" << std::endl;
	model model_pde_sigma(100., sigma, 0.05, 1, 1./10., pp, cond);
	std::cout << "ok5" << std::endl;
	
	solver_edp solver_model(model_pde);
	std::cout << "ok6" << std::endl;
	
	std::vector<double> sol(solver_model.solve_pde());
	
	std::vector<std::vector<double>> mat(3, std::vector<double>(3));
	std::vector<double> d(3);
	
	for(int i=0; i<3; i++)
	{
		d[i] = 3-i;
	}
	
	mat[0][0] = 1;
	mat[1][1] = 2;
	mat[2][2] = 3;
	
	mat[0][1] = -1;
	mat[1][2] = -2;
	
	mat[1][0] = 1;
	mat[2][1] = 2;
	
	std::vector<double> test_inv(solver_model.product_inverse(mat, d));
	
	for(int i=0; i<test_inv.size();i++)
	{
		std::cout << test_inv[i] << ",";
	}
	
	return 0;
}

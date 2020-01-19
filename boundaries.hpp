#ifndef BOUND_HPP
#define BOUND_HPP
#include <vector>
#include "model.hpp"
#include "mesh.hpp"

class bound
{
public:

	bound(payoff f, mesh grille, std::string method, std::vector<std::vector<double>> conditions);
	bound(payoff f, mesh grille, std::string method);
	bound(payoff f, mesh grille, std::vector<std::vector<double>> conditions);
	bound(payoff f, mesh grille);
	
	void get_boundaries(double ri, double ri1, double sigma0, double sigma1, double T, double dt, double j, std::vector<double>& sol);
	
private:
	void get_boundaries_cdt(double ri, double ri1, double sigma0, double sigma1, double T, double dt, double j, std::vector<double>& sol);
	void get_boundaries_nocdt(double ri, double ri1,double sigma0, double sigma1, double T, double	dt, double j, std::vector<double>& sol);
	
	payoff b_f;
	mesh b_mesh;
	
	std::string b_method;
	std::vector<double> b_strikes;
	std::vector<std::vector<double>> b_conditions;
	
};
	

#endif
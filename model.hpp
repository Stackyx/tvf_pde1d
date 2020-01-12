#ifndef MODEL_HPP
#define MODEL_HPP
#include "closed_form.hpp"

class model
{
public:
	model(const int& n, const double& sigma, const double& r, const double& dt, payoff& f);
	
	std::vector<std::vector<double>> pde_matrix(const int& i);
	std::vector<std::vector<double>> pde_matrix_to_inverse(const int& i);
	
private:
	
	double m_sigma;
	double m_r;
	double m_dt;
	int m_n;
	payoff m_f;

};

#endif

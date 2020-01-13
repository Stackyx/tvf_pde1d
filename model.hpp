#ifndef MODEL_HPP
#define MODEL_HPP
#include "closed_form.hpp"

class model
{
public:
	model(const int& n_x, const int& n_t, const double& sigma, const double& r, const double& dt, payoff& f);
	model(const int& n_x, const int& n_t, const std::vector<std::vector<double>>& sigma, const double& r, const double& dt, payoff& f);
	model(const int& n_x, const int& n_t, const double& sigma, const std::vector<double>& r, const double& dt, payoff& f);
	model(const int& n_x, const int& n_t, const std::vector<std::vector<double>>& sigma, const std::vector<double>& r, const double& dt, payoff& f);

	std::vector<std::vector<double>> pde_matrix(const int& i);
	std::vector<std::vector<double>> pde_matrix_to_inverse(const int& i);
	double getS();
	double getS2();
	double get_vol(const int& i, const int& j);
	double get_r(const int&i);
	
private:
	
	std::vector<double> m_cdt;
	double m_initS;
	double m_Smin;
	double m_Smax;
	std::vector<std::vector<double>> m_sigma;
	std::vector<double> m_r
	double m_dx;
	double m_dt;
	int m_n;
	int m_t;
	payoff m_f;

};

#endif

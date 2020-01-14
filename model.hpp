#ifndef MODEL_HPP
#define MODEL_HPP
#include "closed_form.hpp"


class model
{
public:
	model(const double& S0, const double& sigma, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}});
	model(const double& S0, const std::vector<double>& sigma, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}});
	model(const double& S0, const double& sigma, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}});
	model(const double& S0, const std::vector<double>& sigma, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}});

	std::vector<std::vector<double>> pde_matrix(const int& i);
	std::vector<std::vector<double>> pde_matrix_to_inverse(const int& i);
	std::vector<std::vector<double>> getDirichelet();
	
	payoff m_f;
	
	double getS();
	double getS2();
	double get_vol(const int& i);
	double get_r(const int&i);
	
	double m_T;
	double m_Smin;
	double m_dx;
	double m_Smax;
	
	int m_nt;
	int m_nx;
	
private:
	
	std::vector<std::vector<double>> m_cdt;
	double m_initS;
	std::vector<double> m_sigma;
	std::vector<double> m_r;

	double m_dt;
	
	double m_theta;
	
	std::vector<std::vector<double>> get_conditions(std::vector<std::vector<double>> conditions);
	std::vector<double> getStrike();
	std::string getName();
	std::function<double(double)> model::getpayoff();

};

std::vector<double> getRow(std::vector<std::vector<double>> mat, int i);

#endif

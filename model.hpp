#ifndef MODEL_HPP
#define MODEL_HPP
#include "payoff.hpp"


class model
{
public:
	model(const double& S0, const double& sigma, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f);
	model(const double& S0, const std::vector<double>& sigma, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f);
	model(const double& S0, const double& sigma, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f);
	model(const double& S0, const std::vector<double>& sigma, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f);
	
	model(const double& S0, const std::vector<std::vector<double>>& sigma, const double& S_min_mat, const double& S_max_mat, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f);
	model(const double& S0, const std::vector<std::vector<double>>& sigma, const double& S_min_mat, const double& S_max_mat, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f);
	
	double getSmax();
	double getSmin();
	std::vector<double> get_vol_col(const int& i);
	double get_r(const int&i);
	double get_dx();
	std::vector<std::vector<double>> getSigma();
	
	
private:
	
	friend class solver_edp;
	friend class bound;
	
	double m_dt;
	double m_T;
	double m_Smin;
	double m_dx;
	double m_Smax;
	double m_theta;
	double m_initS;

	int m_nt;
	int m_nx;
	
	std::vector<double> m_r;
	std::vector<std::vector<double>> m_sigma;
	
	std::string getName();
	std::function<double(double)> getpayoff();
	std::vector<double> getStrike();
	
	payoff m_f;
	//std::vector<std::vector<double>> resize_sigma(const double& S_min_mat, const double& S_max_mat);
};

std::vector<double> getCol(std::vector<std::vector<double>> mat, int i);
std::vector<double> getRow(std::vector<std::vector<double>> mat, int i);

#endif

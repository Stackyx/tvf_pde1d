#ifndef MODEL_HPP
#define MODEL_HPP
#include "payoff.hpp"


class model
{
public:
	model(double sigma, double r, int n_t, int n_x, double theta, payoff& f);
	model(const std::vector<double>& sigma, double r, int n_t, int n_x, double theta, payoff& f);
	model(double sigma, const std::vector<double>& r, double T, int n_t, int n_x, double theta, payoff& f);
	model(const std::vector<double>& sigma, const std::vector<double>& r, int n_t, int n_x, double theta, payoff& f);
	model(const std::vector<std::vector<double>>& sigma, const double& r, const int& n_t, double theta, payoff& f);
	model(const std::vector<std::vector<double>>& sigma, const std::vector<double>& r, const double theta, payoff& f);
	
	std::vector<double> get_vol_col(const int& i);
	double get_r(const int&i);
	std::vector<std::vector<double>> getSigma();
	
	
private:
	
	friend class solver_edp;
	friend class bound;
	
	double m_theta;
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

#ifndef MODEL_HPP
#define MODEL_HPP
#include "closed_form.hpp"


class model
{
public:
	model(const double& S0, const double& sigma, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}}, std::string method = "Dirichelet");
	model(const double& S0, const std::vector<double>& sigma, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}}, std::string method = "Dirichelet");
	model(const double& S0, const double& sigma, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}}, std::string method = "Dirichelet");
	model(const double& S0, const std::vector<double>& sigma, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}}, std::string method = "Dirichelet");
	
	model(const double& S0, const std::vector<std::vector<double>>& sigma, const double& S_min_mat, const double& S_max_mat, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}}, std::string method = "Dirichelet");
	model(const double& S0, const std::vector<std::vector<double>>& sigma, const double& S_min_mat, const double& S_max_mat, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions = {{0, 0}, {0,0}}, std::string method = "Dirichelet");
	
	std::vector<std::vector<double>> getDirichelet();
	std::vector<std::vector<double>> getNeumann();
	
	double getSmax();
	double getSmin();
	double get_vol(const int& i);
	double get_r(const int&i);
	double get_dx();
	

	
private:
	
	friend class solver_edp;
	
	std::vector<std::vector<double>> m_cdt;
	double m_initS;
	std::vector<std::vector<double>> m_sigma;
	std::vector<double> m_r;

	double m_dt;
	double m_T;
	double m_Smin;
	double m_dx;
	double m_Smax;
	int m_nt;
	int m_nx;
	double m_theta;
	std::string m_method;
	payoff m_f;
	
	std::vector<std::vector<double>> get_conditions(std::vector<std::vector<double>> conditions, std::string method);
	std::vector<double> getStrike();
	std::string getName();
	std::function<double(double)> model::getpayoff();

};

std::vector<double> getRow(std::vector<std::vector<double>> mat, int i);

#endif

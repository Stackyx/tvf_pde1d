#ifndef solver_HPP
#define solver_HPP
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "model.hpp"

model::model(const double& S0, const double& sigma, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions, std::string method)
	: m_nt(n_t), m_nx(n_x), m_T(T), m_f(f), m_initS(S0), m_theta(theta), m_method(method)
	{
		m_dt = T/n_t;
		
		m_sigma.resize(m_nx, std::vector<double>(m_nt));
		
		for (int j = 0; j<m_nt; ++j)
		{
			for(int i=0;i<m_nx;++i)
			{
				m_sigma[i][j] = sigma;
			}
		}
		
		m_r.resize(m_nt);
		
		for(int i=0;i<m_nt;++i)
		{
			m_r[i] = r;
		}
		
		m_Smin = log(m_initS) - 5 * sigma * pow(m_T, 0.5);
		m_Smax = log(m_initS) + 5 * sigma * pow(m_T, 0.5);
		m_dx = (m_Smax - m_Smin)/m_nx;
	
		m_cdt = get_conditions(conditions, method);
	}
	
model::model(const double& S0, const std::vector<double>& sigma, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions, std::string method)
	: m_nt(n_t), m_nx(n_x), m_T(T), m_f(f), m_initS(S0), m_theta(theta)
	{
		
		m_dt = T/n_t;
		double average_vol = accumulate(sigma.begin(), sigma.end(), 0.0)/sigma.size(); 
		//double max_vol = *std::max_element(m_sigma.begin(), m_sigma.end());
		m_Smin = log(m_initS) - 5 * average_vol * pow(m_T, 0.5);
		m_Smax = log(m_initS) + 5 * average_vol * pow(m_T, 0.5);
		m_dx = (m_Smax - m_Smin)/m_nx;
		
		m_r.resize(m_nt);
		
		for(int i=0;i<m_nt;++i)
		{			
			m_r[i] = r;
		}
		
		m_sigma.resize(m_nx,std::vector<double>(m_nt));
	
		for (int i = 0; i<m_nx; ++i)
		{
			std::copy(sigma.begin(), sigma.end(), m_sigma[i].begin());
		}
		
	
		m_cdt = get_conditions(conditions, method);
	}

model::model(const double& S0, const double& sigma, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions, std::string method)
	: m_r(r), m_nt(n_t), m_nx(n_x), m_T(T), m_f(f), m_initS(S0), m_theta(theta)
	{	
		
		m_dt = T/n_t;
		
		m_Smin = log(m_initS) - 5 * sigma * pow(m_T, 0.5);
		m_Smax = log(m_initS) + 5 * sigma * pow(m_T, 0.5);
		m_dx = (m_Smax - m_Smin)/m_nx;
		
		m_sigma.resize(m_nx, std::vector<double>(m_nt));
		
		for (int j = 0; j<m_nt; ++j)
		{
			for(int i=0;i<m_nx;++i)
			{
				m_sigma[i][j] = sigma;
			}
		}
		
		m_cdt = get_conditions(conditions, method);
	}

model::model(const double& S0, const std::vector<double>& sigma, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions, std::string method)
	: m_nt(n_t), m_nx(n_x), m_T(T), m_f(f), m_initS(S0), m_r(r), m_theta(theta)
{
	
	m_dt = T/n_t;
	
	double average_vol = accumulate(sigma.begin(), sigma.end(), 0.0)/sigma.size(); 
	//double max_vol = *std::max_element(m_sigma.begin(), m_sigma.end());
	
	m_Smin = log(m_initS) - 5 * average_vol * pow(m_T, 0.5);
	m_Smax = log(m_initS) + 5 * average_vol * pow(m_T, 0.5);
	m_dx = (m_Smax - m_Smin)/m_nx;
	
	
	m_sigma.resize(m_nx,std::vector<double>(m_nt));
	
	for (int i = 0; i<m_nx; ++i)
	{
		std::copy(sigma.begin(), sigma.end(), m_sigma[i].begin());
	}
	
	m_cdt = get_conditions(conditions, method);
}

model::model(const double& S0, const std::vector<std::vector<double>>& sigma, const double& S_min_mat, const double& S_max_mat, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions, std::string method)
	: m_nt(n_t), m_nx(n_x), m_T(T), m_f(f), m_initS(S0), m_sigma(sigma), m_theta(theta)
{
	
	m_dt = T/n_t;
	m_Smin = log(S_min_mat);
	m_Smax = log(S_max_mat);
	
	size_t size_row_sigma = m_sigma.size();
	m_dx = (S_max_mat - S_min_mat)/size_row_sigma;
	
	m_r.resize(m_nt);
	
	for(int i=0;i<m_nt;++i)
	{			
		m_r[i] = r;
	}
	
	m_cdt = get_conditions(conditions, method);
}

model::model(const double& S0, const std::vector<std::vector<double>>& sigma, const double& S_min_mat, const double& S_max_mat, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions, std::string method)
	: m_nt(n_t), m_nx(n_x), m_T(T), m_f(f), m_initS(S0), m_sigma(sigma), m_r(r), m_theta(theta)
{
	
	m_dt = T/n_t;
	
	
	size_t size_row_sigma = m_sigma.size();
	//double dx_sigma = (S_max_mat - S_min_mat)/size_row_sigma;
	//int i_initS = (m_initS - S_min_mat)/dx_sigma; 
	//double average_vol_t = accumulate(getRow(m_sigma,i_initS).begin(), getRow(m_sigma,i_initS).end(), 0.0)/m_sigma[0].size(); 
	
	m_Smin = log(S_min_mat);
	m_Smax = log(S_max_mat);
	m_dx = (S_max_mat - S_min_mat)/size_row_sigma;
	
	
	m_cdt = get_conditions(conditions, method);
}

// std::vector<std::vector<double>> model::resize_sigma(const double& S_min_mat, const double& S_max_mat)
// {
	// if (m_Smax>S_max_mat)
	// {
		// int n_to_add = (m_Smax - S_max_mat)/m_dx;
		// int size_row = (sizeof(m_sigma)/sizeof(m_sigma[0]));
		// m_sigma.resize(n_to_add + size_row);
		
		// // peut surement etre simplifié avec STL
		// for (int j = 0; j<m_sigma[0].size(); ++j)
		// {
			// for (int i=size_row; i<size_row+n_to_add; ++i)
			// {
				// m_sigma[i][j] = m_sigma[size_row-1][j];
			// }
		// }
	
	// }
	
	
	
// }

double model::getSmax()
{
	return m_Smax;
}
double model::getSmin()
{
	return m_Smin;
}

std::vector<double> model::get_vol_col(const int& i)
{
	return getCol(m_sigma, i);
}

double model::get_r(const int&i)
{
	return m_r[i];
}

double model::get_dx()
{
	return m_dx;
}

std::vector<double> model::getStrike()
{
	return m_f.getparameters();
}

std::string model::getName()
{
	return m_f.getname();
}

std::function<double(double)> model::getpayoff()
{
	return m_f.getpayoff();
}

std::vector<std::vector<double>>  model::getSigma()
{
	return m_sigma;
}


std::vector<double> getCol(std::vector<std::vector<double>> mat, int i)
{
	std::vector<double> temp;
	temp.resize(mat.size());
	
	for (int j = 0; j < (mat.size()); ++j) {
		temp[j] = mat[j][i];
	}
	
	return temp;
}


std::vector<std::vector<double>> model::getDirichelet()
{
	std::vector<double> cdt = getStrike();
	std::vector<std::vector<double>> new_cdt;
	std::vector<double> uppercdt;
	std::vector<double> lowercdt;
	
	new_cdt.resize(m_nt, std::vector<double>(getStrike().size()));
	uppercdt.resize(m_nt);
	lowercdt.resize(m_nt);
	
	for (int j = 0; j<m_nt; ++j)
	{
		for (int i = 0; i<getStrike().size(); ++i)
		{
			new_cdt[j][i] = cdt[i] * exp(- m_r[j] * m_dt*j);

		}
		
		uppercdt[j] = payoff(getName(), new_cdt[j]).getpayoff()(exp(m_Smax));
		lowercdt[j] = payoff(getName(), new_cdt[j]).getpayoff()(exp(m_Smin));
		
	}
	
	return {lowercdt, uppercdt};
	
}


std::vector<std::vector<double>> model::getNeumann()
{
	//a modifié
	std::vector<double> cdt = getStrike();
	std::vector<std::vector<double>> new_cdt;
	std::vector<double> uppercdt;
	std::vector<double> lowercdt;
	
	new_cdt.resize(m_nt, std::vector<double>(getStrike().size()));
	uppercdt.resize(m_nt);
	lowercdt.resize(m_nt);
	
	for (int j = 0; j<m_nt; ++j)
	{
		for (int i = 0; i<getStrike().size(); ++i)
		{
			new_cdt[j][i] = cdt[i] * exp(- m_r[j] * m_dt*j);

		}
		
		uppercdt[j] = payoff(getName(), new_cdt[j]).getpayoff()(exp(m_Smax));
		lowercdt[j] = payoff(getName(), new_cdt[j]).getpayoff()(exp(m_Smin));
		
	}
	
	return {lowercdt, uppercdt};
	
}


std::vector<std::vector<double>> model::get_conditions(std::vector<std::vector<double>> conditions, std::string method)
{	
	
	std::vector<std::vector<double>> c = {{0, 0}, {0,0}};
	if (std::equal(conditions[0].begin(), conditions[0].end(), c[0].begin()) && std::equal(conditions[1].begin(), conditions[1].end(), c[1].begin()))
	{
		if (CaseSensitiveIsEqual(method, "Dirichelet"))
		{
			std::vector<std::vector<double>> drchlt = getDirichelet();
			conditions.resize(2, std::vector<double>(m_nt));
			conditions = drchlt;
		}
		else if (CaseSensitiveIsEqual(method, "Newmann"))
		{
			std::vector<std::vector<double>> Nwn = getNeumann();
			conditions.resize(2, std::vector<double>(m_nt));
			conditions = Nwn;
		}
		else
		{
			std::cout<< "Issue with the name of the method" << std::endl;
		}
	}
	
	return conditions;
}


#endif

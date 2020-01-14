#ifndef solver_HPP
#define solver_HPP
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "model.hpp"

model::model(const double& S0, const double& sigma, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions)
	: m_nt(n_t), m_nx(n_x), m_T(T), m_f(f), m_initS(S0), m_theta(theta)
	{
		m_dt = T/n_t;
		
		m_sigma.resize(m_nt);
		
		for(int i=0;i<m_nt;++i)
		{
			m_sigma[i] = sigma;
		}
		
		m_r.resize(m_nt);
		
		for(int i=0;i<m_nt;++i)
		{
			m_r[i] = r;
		}
		
		m_cdt = get_conditions(conditions);
	}
	
model::model(const double& S0, const std::vector<double>& sigma, const double& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions)
	: m_nt(n_t), m_nx(n_x), m_T(T), m_f(f), m_initS(S0), m_sigma(sigma), m_theta(theta)
	{
		
		m_dt = T/n_t;

		m_r.resize(m_nt);
		
		for(int i=0;i<m_nt;++i)
		{			
			m_r[i] = r;
		}
		
		m_cdt = get_conditions(conditions);
	}

model::model(const double& S0, const double& sigma, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions)
	: m_r(r), m_nt(n_t), m_nx(n_x), m_T(T), m_f(f), m_initS(S0), m_theta(theta)
	{	
		
		m_dt = T/n_t;
		
		m_sigma.resize(m_nt);
		
		for(int i=0;i<m_nt;++i)
		{
			m_sigma[i] = sigma;
		}
		
		m_cdt = get_conditions(conditions);
	}

model::model(const double& S0, const std::vector<double>& sigma, const std::vector<double>& r, const double& T, const int& n_t, const int& n_x, const double& theta, payoff& f, std::vector<std::vector<double>> conditions)
	: m_nt(n_t), m_nx(n_x), m_T(T), m_f(f), m_initS(S0), m_sigma(sigma), m_r(r), m_theta(theta)
{
	
	m_dt = T/n_t;
		
	m_cdt = get_conditions(conditions);
}

double model::getS()
{
	return m_Smax;
}
double model::getS2()
{
	return m_Smin;
}

double model::get_vol(const int& i)
{
	return m_sigma[i];
}

double model::get_r(const int&i)
{
	return m_r[i];
}

std::vector<std::vector<double>> model::pde_matrix(const int& i)
{	
	std::vector<std::vector<double>> mat(3, std::vector<double>(m_nx));
	
	mat[1][0] = 1;
	mat[1][m_nx-1] = 1;
	
	double sigma(get_vol(i));
	double r(get_r(i));
	
	for(int j = 1; j<m_nx-1; ++j)
	{
		mat[1][j] = 1 - (1-m_theta)*(pow(sigma/m_dx,2) + r);
		mat[0][j] = (1-m_theta)/(2*m_dx)*(pow(sigma,2)/m_dx + pow(sigma,2)/2 - r);
		mat[2][j] = (1-m_theta)/(2*m_dx)*(pow(sigma,2)/m_dx - pow(sigma,2)/2 + r);
	}
	
	return mat;
}

std::vector<std::vector<double>> model::pde_matrix_to_inverse(const int& i)
{
	std::vector<std::vector<double>> mat(3, std::vector<double>(m_nx));
	
	// Conditions initiales et terminales sur X :
	
	mat[1][0] = 1;
	mat[1][m_nx-1] = 1;
	
	double sigma(get_vol(i));
	double r(get_r(i));
	
	for(int j = 1; j<m_nx-1; ++j)
	{		
		mat[1][j] = 1+m_dt*m_theta*(pow(sigma/m_dx,2) + r);
		mat[0][j] = m_dt*m_theta/(2*m_dx)*(-pow(sigma,2)/m_dx - pow(sigma,2)/2 + r);
		mat[2][j] = m_dt*m_theta/(2*m_dx)*(-pow(sigma,2)/m_dx + pow(sigma,2)/2 - r);
	}
	
	return mat;
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

std::vector<double> getRow(std::vector<std::vector<double>> mat, int i)
{
	std::vector<double> temp;
	temp.resize(mat[i].size());
	
	for (int j = 0; j < mat[i].size(); ++j) {
		temp[j] = mat[i][j];
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
			new_cdt[j][i] = cdt[i] * exp(- m_r[j] * (m_dt*j));

		}
		
		uppercdt[j] = payoff(getName(), getRow(new_cdt, j)).getpayoff()(exp(m_Smax));
		lowercdt[j] = payoff(getName(), getRow(new_cdt, j)).getpayoff()(exp(m_Smin));
		
	}
	
	return {lowercdt, uppercdt};
	
}


std::vector<std::vector<double>> model::get_conditions(std::vector<std::vector<double>> conditions)
{
	
	double average_vol = accumulate(m_sigma.begin(), m_sigma.end(), 0.0)/m_sigma.size(); 
	double max_vol = *std::max_element(m_sigma.begin(), m_sigma.end());
	
	m_Smin = log(m_initS) - 5 * average_vol * pow(m_T, 0.5);
	m_Smax = log(m_initS) + 5 * average_vol * pow(m_T, 0.5);
	
	m_dx = (m_Smax - m_Smin)/m_nx;
	
	
	std::vector<std::vector<double>> c = {{0, 0}, {0,0}};
	if (std::equal(conditions[0].begin(), conditions[0].end(), c[0].begin()) && std::equal(conditions[1].begin(), conditions[1].end(), c[1].begin()))
	{
		std::vector<std::vector<double>> drchlt = getDirichelet();
		conditions.resize(2, std::vector<double>(m_nt));
		conditions = drchlt;
	}
	
	return conditions;
}


#endif

#ifndef solver_HPP
#define solver_HPP
#include <vector>
#include <iostream>
#include "model.hpp"

model::model(const double& S0, const double& sigma, const double& r, const double& T, const double& dt, payoff& f, std::vector<double> conditions)
	: m_dt(dt), m_T(T), m_f(f), m_initS(S0)
	{
		m_nt = std::round(T/dt);
		
		for(int i=0;i<m_nx;++i)
		{
			for(int j=0;j<m_nx;++j)
			{
				m_sigma[i][j] = sigma;
			}
		}
		
		m_r.resize(m_nt);
		
		for(int i=0;i<m_nt;++i)
		{
			m_r[i] = r;
		}
		
		m_cdt = get_conditions(conditions);
	}
	
model::model(const double& S0, const std::vector<std::vector<double>>& sigma, const double& r, const double& T, const double& dt, payoff& f, std::vector<double> conditions)
	: m_dt(dt), m_T(T), m_f(f), m_initS(S0), m_sigma(sigma)
	{
		m_nt = std::round(T/dt);
		
		m_r.resize(m_nt);
		
		for(int i=0;i<m_nt;++i)
		{			
			m_r[i] = r;
		}
		
		m_cdt = get_conditions(conditions);
	}

model::model(const double& S0, const double& sigma, const std::vector<double>& r, const double& T, const double& dt, payoff& f, std::vector<double> conditions)
	: m_r(r), m_dt(dt), m_T(T), m_f(f), m_initS(S0)
	{	
		
		m_nt = std::round(T/dt);
		
		std::vector<std::vector<double>> m_sigma(m_nx, std::vector<double>(m_nx));
		
		for(int i=0;i<m_nx;++i)
		{
			for(int j=0;j<m_nx;++j)
			{
				m_sigma[i][j] = sigma;
			}
		}
		
		m_cdt = get_conditions(conditions);
	}

model::model(const double& S0, const std::vector<std::vector<double>>& sigma, const std::vector<double>& r, const double& T, const double& dt, payoff& f, std::vector<double> conditions)
	: m_dt(dt), m_T(T), m_f(f), m_initS(S0), m_sigma(sigma), m_r(r)
{
	
	m_nt = std::round(T/dt);
		
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

double model::get_vol(const int& i, const int& j)
{
	return m_sigma[i][j];
}

double model::get_r(const int&i)
{
	return m_r[i];
}

std::vector<std::vector<double>> model::pde_matrix(const int& i)
{	
	std::vector<std::vector<double>> mat(3, std::vector<double>(m_nx));
	
	mat[1][0] = 1;
	mat[0][0] = 0;
	mat[2][m_nx-2] = 0;
	
	for(int j = 1; j<m_nx-2; ++j)
	{
		mat[1][j] = 0;
		mat[0][j] = 0;
		mat[2][j] = 0;
	}
	
	return mat;
}

std::vector<std::vector<double>> model::pde_matrix_to_inverse(const int& i)
{
	std::vector<std::vector<double>> mat(3, std::vector<double>(m_nx));
	
	mat[1][0] = 1;
	mat[0][0] = 0;
	mat[2][m_nx-2] = 0;
	
	for(int j = 1; j<m_nx-2; ++j)
	{
		mat[1][j] = 0;
		mat[0][j] = 0;
		mat[2][j] = 0;
	}
	
	return mat;
}

std::vector<double> model::get_conditions(std::vector<double> conditions)
{
	
	m_Smin = exp(log(m_initS) - 5 * get_vol(m_nt-1,0) * pow(m_T, 0.5));
	m_Smax = exp(log(m_initS) + 5 * get_vol(m_nt-1, m_nx-1) * pow(m_T, 0.5));
	
	m_dx = (m_Smax - m_Smin) * m_dt;
	m_nx = int ((m_Smax - m_Smin)/m_dx);
	
	std::vector<double> c = { 0, 9999999 };
	if (std::equal(conditions.begin(), conditions.end(), c.begin()))
	{
		conditions[0] = m_f.getpayoff()(m_Smin);
		conditions[1] = m_f.getpayoff()(m_Smax);
	}
	
	return conditions;
}

#endif

#ifndef solver_HPP
#define solver_HPP
#include <vector>
#include <iostream>
#include "model.hpp"

model::model(const double& S0, const int& n, const double& sigma, const double& r, const double& T, payoff& f, std::vector<double> conditions)
	: m_sigma(sigma), m_r(r), m_n(n), m_f(f), m_initS(S0)
	{
		m_dt = 1.0 / 252.0;
		m_Smin = exp(log(m_initS) - 5 * m_sigma * pow(T, 0.5));
		m_Smax = exp(log(m_initS) + 5 * m_sigma * pow(T, 0.5));
		m_dx = (m_Smax - m_Smin) * m_dt;
		
		std::vector<double> c = { 0, 9999999 };
		if (std::equal(conditions.begin(), conditions.end(), c.begin()))
		{
			conditions[0] = f.getpayoff()(m_Smin);
			conditions[1] = f.getpayoff()(m_Smax);
		}
		
		m_r.resize(m_t);
		
		for(int i=0;i<m_t;++i)
		{
			m_r[i] = r;
		}
	}

		m_cdt = conditions;

model::model(const int& n_x, const int& n_t, const std::vector<std::vector<double>>& sigma, const double& r, const double& dt, payoff& f)
	: m_n(n_x), m_dt(dt), m_f(f), m_t(n_t)
	{
		m_r.resize(m_t);
		
		for(int i=0;i<m_t;++i)
		{			
			m_r[i] = r;
		}
	}

model::model(const int& n_x, const int& n_t, const double& sigma, const std::vector<double>& r, const double& dt, payoff& f)
	: m_r(r), m_n(n_x), m_dt(dt), m_f(f), m_t(n_t)
	{
		std::vector<std::vector<double>> m_sigma(m_n, std::vector<double>(m_n));
		
		for(int i=0;i<m_n;++i)
		{
			for(int j=0;j<m_n;++j)
			{
				m_sigma[i][j] = sigma;
			}
		}
	}

double model::getS()
{
	return m_Smax;
}
double model::getS2()
{
	return m_Smin;
}

std::vector<std::vector<double>> model::pde_matrix(const int& i)
{	
	std::vector<std::vector<double>> mat(m_n, std::vector<double>(m_n));
	
	
	
	for(int j = 1; j<m_n-1; ++j)
	{
		mat[j][j] = 1-1./2.*get_r(i)*i-1./2.*pow(get_vol(i,j)*i,2);
		mat[j][j-1] = 1./4.*pow(get_vol(i,j)*i,2);
		mat[j][j+1] = 1./2.*get_r(i)*i+1./4.*pow(get_vol(i,j)*i,2);
	}
	
	return mat;
}

std::vector<std::vector<double>> model::pde_matrix_to_inverse(const int& i)
{
	std::vector<std::vector<double>> mat(m_n, std::vector<double>(m_n));

	for(int j = 1; j<m_n-1; ++j)
	{
		mat[j][j] = m_dt*(get_r(i) - 1./m_dt + 1./2.*get_r(i)*i+1./2.*pow(get_vol(i,j)*i, 2));
		mat[j][j-1] = -1./4.*m_dt*pow(get_vol(i,j)*i, 2);
		mat[j][j+1] = m_dt*(-1./2.*get_r(i)*i-1./4.*pow(get_vol(i,j)*i,2));
	}
	
	return mat;
}

#endif

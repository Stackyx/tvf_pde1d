#ifndef solver_HPP
#define solver_HPP
#include <vector>
#include <iostream>
#include "model.hpp"

model::model(const int& n_x, const int& n_t, const double& sigma, const double& r, const double& dt, payoff& f)
	: m_n(n_x), m_dt(dt), m_f(f), m_t(n_t)
	{
		
		std::vector<std::vector<double>> m_sigma(m_n, std::vector<double>(m_n));
		
		for(int i=0;i<m_n;++i)
		{
			for(int j=0;j<m_n;++j)
			{
				m_sigma[i][j] = sigma;
			}
		}
		
		m_r.resize(m_t);
		
		for(int i=0;i<m_t;++i)
		{
			m_r[i] = r;
		}
	}

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

model::model(const int& n_x, const int& n_t, const std::vector<std::vector<double>>& sigma, const std::vector<double>& r, const double& dt, payoff& f)
	: m_sigma(sigma), m_r(r), m_n(n_x), m_dt(dt), m_f(f), m_t(n_t)
	{
	}

double model::get_vol(const int& i, const int& j)
{
	return m_sigma[i][j];
}

double model::get_r(const int& i)
{
	return m_r[i];
}

std::vector<std::vector<double>> model::pde_matrix(const int& i)
{	
	std::vector<std::vector<double>> mat(m_n, std::vector<double>(m_n));
	
	// Need to apply boundaries conditions
	
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
	
	// Need to apply boundaries conditions
	
	for(int j = 1; j<m_n-1; ++j)
	{
		mat[j][j] = m_dt*(get_r(i) - 1./m_dt + 1./2.*get_r(i)*i+1./2.*pow(get_vol(i,j)*i, 2));
		mat[j][j-1] = -1./4.*m_dt*pow(get_vol(i,j)*i, 2);
		mat[j][j+1] = m_dt*(-1./2.*get_r(i)*i-1./4.*pow(get_vol(i,j)*i,2));
	}
	
	return mat;
}

#endif
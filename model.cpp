#ifndef solver_HPP
#define solver_HPP
#include <vector>
#include <iostream>
#include "model.hpp"

model::model(const int& n, const double& sigma, const double& r, const double& dt, payoff& f)
	: m_sigma(sigma), m_r(r), m_n(n), m_dt(dt), m_f(f)
	{
	}
	
std::vector<std::vector<double>> model::pde_matrix(const int& i)
{	
	std::vector<std::vector<double>> mat(m_n, std::vector<double>(m_n));
	
	// Need to apply boundaries conditions
	
	for(int j = 1; j<m_n-1; ++j)
	{
		mat[j][j] = 1-1./2.*m_r*i-1./2.*pow(m_sigma*i,2);
		mat[j][j-1] = 1./4.*pow(m_sigma*i,2);
		mat[j][j+1] = 1./2.*m_r*i+1./4.*pow(m_sigma*i,2);
	}
	
	return mat;
}

std::vector<std::vector<double>> model::pde_matrix_to_inverse(const int& i)
{
	std::vector<std::vector<double>> mat(m_n, std::vector<double>(m_n));
	
	// Need to apply boundaries conditions
	
	for(int j = 1; j<m_n-1; ++j)
	{
		mat[j][j] = m_dt*(m_r - 1./m_dt + 1./2.*m_r*i+1./2.*pow(m_sigma*i, 2));
		mat[j][j-1] = -1./4.*m_dt*pow(m_sigma*i, 2);
		mat[j][j+1] = m_dt*(-1./2.*m_r*i-1./4.*pow(m_sigma*i,2));
	}
	
	return mat;
}

#endif
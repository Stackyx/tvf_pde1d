#ifndef SOLVER_HPP
#define SOLVER_HPP
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "model.hpp"

model::model(double sigma, double r, int n_t, int n_x)
	{	
		
		if (n_x % 2 == 0)
		{
			n_x += 1;
		}
		
		m_sigma.resize(n_x, std::vector<double>(n_t));
		m_r.resize(n_t);
		
		for (int j = 0; j<n_t; ++j)
		{
			for(int i=0;i<n_x;++i)
			{
				m_sigma[i][j] = sigma;
			}
			
			m_r[j] = r;
		}
		
	}
	
model::model(const std::vector<double>& sigma, double r, int n_t, int n_x)
	{
		
		if (n_x % 2 == 0)
		{
			n_x += 1;
		}
		
		m_r.resize(n_t);
		
		for(int i=0;i<n_t;++i)
		{			
			m_r[i] = r;
		}
		
		m_sigma.resize(n_x,std::vector<double>(n_t));
	
		for (int i = 0; i<n_x; ++i)
		{
			std::copy(sigma.begin(), sigma.end(), m_sigma[i].begin());
		}

	}

model::model(double sigma, const std::vector<double>& r, double T, int n_t, int n_x)
	: m_r(r)
	{	
	
		if (n_x % 2 == 0)
		{
			n_x += 1;
		}
		
		m_sigma.resize(n_x, std::vector<double>(n_t));
		
		for (int j = 0; j<n_t; ++j)
		{
			for(int i=0;i<n_x;++i)
			{
				m_sigma[i][j] = sigma;
			}
		}

	}

model::model(const std::vector<double>& sigma, const std::vector<double>& r, int n_t, int n_x)
	: m_r(r)
{
	
	if (n_x % 2 == 0)
	{
		n_x += 1;
	}
	
	m_sigma.resize(n_x,std::vector<double>(n_t));
	
	for (int i = 0; i<n_x; ++i)
	{
		std::copy(sigma.begin(), sigma.end(), m_sigma[i].begin());
	}

}

model::model(const std::vector<std::vector<double>>& sigma, const double& r, const int& n_t)
	: m_sigma(sigma)
{
	m_r.resize(n_t);
	
	for(int i=0;i<n_t;++i)
	{			
		m_r[i] = r;
	}

}

model::model(const std::vector<std::vector<double>>& sigma, const std::vector<double>& r)
	: m_sigma(sigma), m_r(r)
{
}

std::vector<double> model::get_vol_col(const int& i)
{
	return getCol(m_sigma, i);
}

double model::get_r(const int&i)
{
	return m_r[i];
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

#endif

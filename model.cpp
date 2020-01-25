#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>

#include "tools.hpp"
#include "model.hpp"

model::model(double sigma, double r, int n_t, int n_x)
{	
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
	
	size_check(m_sigma, m_r, n_t, n_x);
		
}
	
model::model(const std::vector<double>& sigma, double r, int n_t, int n_x)
{
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
	
	size_check(m_sigma, m_r, n_t, n_x);

}

model::model(double sigma, const std::vector<double>& r, int n_t, int n_x)
	: m_r(r)
{	
		
	m_sigma.resize(n_x, std::vector<double>(n_t));
	
	for (int j = 0; j<n_t; ++j)
	{
		for(int i=0;i<n_x;++i)
		{
			m_sigma[i][j] = sigma;
		}
	}
		
	size_check(m_sigma, m_r, n_t, n_x);

}

model::model(const std::vector<double>& sigma, const std::vector<double>& r, int n_t, int n_x)
	: m_r(r)
{
	
	m_sigma.resize(n_x,std::vector<double>(n_t));
	
	for (int i = 0; i<n_x; ++i)
	{
		std::copy(sigma.begin(), sigma.end(), m_sigma[i].begin());
	}
	
	size_check(m_sigma, m_r, n_t, n_x);
}

model::model(const std::vector<std::vector<double>>& sigma, double r, int n_t, int n_x)
	: m_sigma(sigma)
{
	m_r.resize(n_t);
	
	for(int i=0;i<n_t;++i)
	{			
		m_r[i] = r;
	}
	
	size_check(m_sigma, m_r, n_t, n_x);
}

model::model(const std::vector<std::vector<double>>& sigma, const std::vector<double>& r, int n_t, int n_x)
	: m_sigma(sigma), m_r(r)
{
	
	size_check(m_sigma, m_r, n_t, n_x);
	
}

void model::get_vol_col(std::vector<double>& mat, int i) const
{	
	for (int j = 0; j < mat.size(); ++j)
	{
		mat[j] = m_sigma[j][i];
	}
}

double model::get_r(int i) const
{
	return m_r[i];
}

double model::get_r_avg() const
{
	return std::accumulate(m_r.begin(), m_r.end(), 0.0)/m_r.size(); // Works ok if dt is constant (it is)
}

std::vector<double> model::get_r() const
{
	return m_r;
}

std::vector<std::vector<double>>  model::getSigma() const
{
	return m_sigma;
}

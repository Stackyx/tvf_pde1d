#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "mesh.hpp"

mesh::mesh(const double& S, const double& T, const int& n_x, const int& n_t, const double& sigma)
	: m_nx(n_x), m_nt(n_t), m_S(S), m_sigma(sigma)
{
	if (m_nx % 2 == 0)
	{
		m_nx += 1;
	}
	
	m_dt = T/(m_nt - 1);
	
	m_Smax = std::log(S) + 5 * sigma * std::sqrt(T);
	m_Smin = std::log(S) - 5 * sigma * std::sqrt(T);
	
	m_dx = (m_Smax - m_Smin)/(m_nx - 1);
}

void mesh::print_mesh()
{
	for(int j=0; j < m_nx; ++j)
	{
		std::cout << std::exp(m_Smin + j*m_dx) << ",";
	}
	
	std::cout << std::endl;
	
	for(int i=0; i < m_nt; ++i)
	{
		std::cout << i*m_dt << std::endl;
	}
}

double mesh::get_dx()
{
	return m_dx;
}

double mesh::get_dt()
{
	return m_dt;
}

double mesh::get_Smax()
{
	return m_Smax;
}

double mesh::get_Smin()
{
	return m_Smin;
}

int mesh::get_nx()
{
	return m_nx;
}

int mesh::get_nt()
{
	return m_nt;
}

double mesh::get_S()
{
	return m_S;
}

double mesh::get_sigma()
{
	return m_sigma;
}

double mesh::get_mat()
{
	return (m_nt-1)*m_dt;
}
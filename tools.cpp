#include "tools.hpp"
#include <exception>
#include <iostream>
#include <exception>

double sqr(double x)
{
	return x*x;
}

void size_check(const std::vector<std::vector<double>>& m_sigma, const std::vector<double>& m_r, int nt, int nx)
{
	std::string msg = "";
	
	if (m_sigma.size() != nx)
	{
		msg = msg + "Sigma size mispecified in x. ";
	}
	if (m_sigma[0].size() != nt)
	{
		msg = msg + "Sigma size mispecified in t. ";
	}
	if (m_r.size() != nt)
	{
		msg = msg + "Rate size mispecified in t.";
	}
	
	if (msg != "")
	{
		throw std::runtime_error(msg);
	}
}
#ifndef UVECTOR_HPP
#define UVECTOR_HPP
#include <vector>
#include "solver.hpp"

solver_edp::solver_edp(const double& stepsT, const double& stepsX, const std::vector<double>& u0, std::vector<std::vector<double>> vol)
	: s_stepsT(stepsT), s_stepsX(stepsX), s_uo(u0), s_vol(vol)
	{
	}

std::vector<std::vector<double>> solver_edp::get_matrix()
{
	std::vector<std::vector<double>> mat(v_vol.size(), 0);
	
	std::vector<std::vector<double>> kappa = get_kappa(const int& j);
	
	for(int i = 0; i < s_step; ++i)
	{
		if (i == 0)
		{
			mat[i][0] = -2/(s_step*s_step);
			mat[i][1] = 1/(s_step*s_step);
		}
		else if (i < s_step-1)
		{
			mat[i][i] = -2/(s_step*s_step);
			mat[i][i-1] = 1/(s_step*s_step);
			mat[i][i+1] = 1/(s_step*s_step);
		}
		else
		{
			mat[i][i] = -2/(s_step*s_step);
			mat[i][i-1] = 1/(s_step*s_step);
		}
	}
	
	return mat;
}

std::vector<std::vector<double>> solver_edp::get_kappa(const int& j)
{
	std::vector<std::vector<double>> kappa(v_vol.size(), 0);
	
	for (int i = 0; i < v_vol.size(); ++i)
	{
		kappa[i][j] = s_stepsT/pow(s_stepsX, 2)*v_vol[i][j]/2
	}
				
	return kappa;

vol_model::vol_model(const std::vector<std::vector<double>>& vol)
	: v_vol(vol)
	{
	}
}

#endif
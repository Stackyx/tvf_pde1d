#include <vector>
#include "solver.hpp"
#include <iostream>

solver_edp::solver_edp(model pde_model)
	: s_pde_model(pde_model)
{
}

std::vector<double> solver_edp::solve_pde()
{	
	std::vector<double> a;
	return a;
}

std::vector<double> solver_edp::product_inverse(std::vector<std::vector<double>> mat, std::vector<double> d)
{
	double w;
	std::vector<double> a(mat[0].size());
	std::vector<double> b(mat[0].size());
	std::vector<double> c(mat[0].size());
	
	std::vector<double> x(mat[0].size());

	for(int i=0; i<b.size(); ++i)
	{
		b[i] = mat[i][i];
		
		if (i==0)
		{
			a[i] = 0;
			c[i] = mat[i][i+1];
		}
		else if(i==b.size()-1)
		{
			c[i] = 0;
			a[i] = mat[i][i-1];
		}
		else
		{
			a[i] = mat[i][i-1];
			c[i] = mat[i][i+1];
		}
	}
	
	for(int i=1; i< mat[0].size(); ++i)
	{
		w = a[i]/b[i-1];
		b[i] = b[i] - w*c[i-1];
		d[i] = d[i] - w*d[i-1];
	}
	
	x[mat[0].size()-1] = d[mat[0].size()-1]/b[mat[0].size()-1];

	for(int i=mat[0].size()-2; i >= 0; --i)
	{
		x[i] = (d[i]-c[i]*x[i+1])/b[i];
	}
	
	return x;
}
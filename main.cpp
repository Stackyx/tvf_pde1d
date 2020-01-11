#include <iostream>
#include "closed_form.hpp"
#include "solver.hpp"

int main(int argc, char* argv[])
{
	
	solver_edp test(1, 10);
	
	std::vector<std::vector<float>> test_mat;
	
	test_mat = test.get_matrix();
	
	for (int i = 0; i < test_mat.size(); ++i)
	{
		for (int j = 0; j < test_mat.size(); ++j)
		{
			std::cout << test_mat[i][j] << " ";
		}
		std::cout << std::endl;
	}
	
	return 0;
}

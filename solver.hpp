#include <vector>

class solver_edp
{
public:
	
	solver_edp(const double& maturity, const int& step, const vol_model& vol);
	std::vector<std::vector<double>> get_matrix();
	std::vector<std::vector<double>> get_kappa();
	
private:
	double s_stepsT;
	double s_stepsX;
	std::vector<std::vector<double>> s_vol;
	std::vector<double> u0;
};

class vol_model
{
	vol_model(const std::vector<std::vector<double>>& vol);
	
private:
	std::vector<std::vector<double>> v_vol;
}
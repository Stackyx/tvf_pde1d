#ifndef PAYOFF_HPP
#define PAYOFF_HPP

#include <vector>
#include <string>
#include <functional>

class payoff
{
public:
	explicit payoff(const std::string& name = "", const std::vector<double>& parameters = { 0 }, std::function<double(double)> fct = [](double d1) { return d1; });
	std::string getname();
	std::function<double(double)> getpayoff();
	std::vector<double>& getparameters();

private:
	std::string m_name;
	std::vector<double> param;
	std::function<double(double)> payoff_fct;
};

bool CaseSensitiveIsEqual(const std::string& str1, const std::string str2);


#endif

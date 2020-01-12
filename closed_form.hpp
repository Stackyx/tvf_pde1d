#ifndef CLOSED_FORM_HPP
#define CLOSED_FORM_HPP

#include <vector>
#include <string>
#include <functional>


double vanilla_payoff(double fwd, double strike, bool is_call);
double bs_time_value(double fwd, double strike, double volatility, double maturity);
double bs_price(double fwd, double strike, double volatility, double maturity, bool is_call);

std::vector<double> vanilla_payoff(const std::vector<double>& fwd, double strike, bool is_call);
std::vector<double> bs_time_value(const std::vector<double>& fwd, double strike, double volatility, double maturity);
std::vector<double> bs_price(const std::vector<double>& fwd, double strike, double volatility, double maturity, bool is_call);

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

bool CaseSensitiveIsEqual(std::string& str1, std::string str2);


#endif

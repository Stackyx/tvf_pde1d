#include "payoff.hpp"

#include <cmath>
#include <limits>
#include <algorithm>
#include <string>
#include <cctype>

payoff::payoff(const std::string& name, const std::vector<double>& parameters, std::function<double(double)> fct)
		:m_name(name), param(parameters)
	{
		std::string s;

		if (CaseSensitiveIsEqual(m_name, "Call")) {
			payoff_fct = [&](double d1) { return std::max(d1 - param[0], 0.0); };
		}
		else if (CaseSensitiveIsEqual(m_name, "Put")) {
			payoff_fct = [&](double d1) { return std::max(param[0] - d1, 0.0); };
		}
		else if (CaseSensitiveIsEqual(m_name, "Straddle")) {
			payoff_fct = [&](double d1) { return std::max(param[0] - d1, d1 - param[0]); };
		}
		else if (CaseSensitiveIsEqual(m_name, "Call Spread")) {
			payoff_fct = [&](double d1) { return std::max(0.0, std::min(d1 - param[0], param[1] - param[0])); };
		}
		else if (CaseSensitiveIsEqual(m_name, "Put Spread")) {
			payoff_fct = [&](double d1) { return std::min(param[1] - param[0], std::max(param[1] - d1, 0.0)); };
		}
		else if (CaseSensitiveIsEqual(m_name, "Strangle")) {
			payoff_fct = [&](double d1) { return std::max(param[0] - d1, 0.0) + std::max(d1 - param[1], 0.0); };
		}
		else {
			payoff_fct = fct;
		}

	}

std::vector<double>& payoff::getparameters()
{
	return param;
}

std::string payoff::getname()
{
	return m_name;
}

std::function<double(double)> payoff::getpayoff()
{
	return payoff_fct;
}

bool CaseSensitiveIsEqual(const std::string& str1, const std::string str2)
{
	return ((str1.size() == str2.size()) && std::equal(str1.begin(), str1.end(), str2.begin(), [](char c1, char c2) {
		return (c1 == c2 || std::toupper(c1) == std::toupper(c2));
		}));
}

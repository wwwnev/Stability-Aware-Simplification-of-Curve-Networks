#pragma once

#include <Eigen/Core>
#include "curve_network.h"


class cn_newton_solver
{
public:
	curve_network m_CN_initial;
	curve_network m_CN_current;

	void init();
	double compute_energy_and_gradient(const Eigen::VectorXd& x, Eigen::VectorXd& grad, bool compute_gradient = true);
	double compute_energy(const Eigen::VectorXd& x);
	void solve(const std::string& s = "", int max_i = 0, double epsilon = -1);
	
	cn_newton_solver(curve_network cn);
};

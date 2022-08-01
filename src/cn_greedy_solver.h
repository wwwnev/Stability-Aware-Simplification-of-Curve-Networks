#pragma once
#include <Eigen/Core>
#include "curve_network.h"

class cn_greedy_solver {
	public:
		curve_network m_cn;
		double m_lambda;

		cn_greedy_solver(const curve_network& cn);
		
		void solve_with_relax(Eigen::VectorXd& alphas, const std::string& outname = "", double budget = 0, std::vector<int> keep = {}, std::vector<int> kick = {});

		void solve_with_load(Eigen::VectorXd& alphas, const std::string& outname = "", double budget = 0, std::vector<int> keep = {}, std::vector<int> kick = {});
};
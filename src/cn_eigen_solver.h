#pragma once
#include <vector>
#include <iostream>
#include <algorithm>
#include "set"
#include <math.h>
#include <numeric>
#include <execution>
#include <chrono>
#include <string>
#include <iostream>
#include <future>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <polyscope/curve_network.h>

#include "constants.h"
#include "aux_tools.h"
#include "curve.h"
#include "curve_network.h"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class cn_eigen_solver {
public:
	curve_network m_cn;
	int m_n;
	int m_nb_curves;
	double m_budget;

	double m_lambda;
	Eigen::VectorXd m_all_lambdas;

	std::vector<SpMat> m_qqt;
	//std::vector<std::vector<SpMat>> m_qqt_intersection;
	std::vector<SpMat> m_qqt_intersection;
	double m_qqt_weight;
	double m_qqt_weight_coefficient;

	std::vector<SpMat> m_identity_reg;

	// this probably obsolete now
	std::map<int, std::pair<int, int>> m_betas_map; // key = beta idx, value = idx of the 2 curves intersecting

	cn_eigen_solver(curve_network& cn, double qqt_weight_coeff = 1.3); // constructor

	void compute_qqt();
	void compute_qqt_intersections();
	void compute_betas_map();
	void compute_identity_regularizers();
	void layer_hessians_fast(const Eigen::VectorXd& alphas, SpMat& H, double qqt_w = 1.0, double w = 1.0);
	void layer_hessians(const Eigen::VectorXd& alphas, SpMat& H, double qqt_w = 1.0, double w = 1.0); // w is stiffness cheat for floating segments.
	void brute_force_solve(Eigen::VectorXd& alphas, const std::string& outname = "");
	void brute_force_solve_nb_curves_budget(Eigen::VectorXd& alphas, int budget = 1, const std::string& outname = "");
	void brute_force_solve_slow(Eigen::VectorXd& alphas, const std::string& outname = "");
	std::vector<double> brute_force_solve_v2(double greedy_lambda, double budget);

	void generate_legal_combination(double budget, const std::vector<int>& kick, const std::vector<int>& keep, Eigen::VectorXd& c);

	void parallel_for_brute_force(std::future<void> futureObj, Eigen::VectorXd& all_lambdas, const Eigen::MatrixXd& all_combs, int greedy_budget, const std::vector<int>& greedy_kick, const std::vector<int>& keep);

	std::vector<double> random_sample_brute_force_v2(double greedy_lambda, double greedy_time, double greedy_budget, Eigen::VectorXd& max_alpha, std::vector<int> greedy_kick, const std::string& outname, std::vector<int> keep);

	void greedy_solve_one_iter(double budget, double& cur_lambda, int& removed_curve, std::vector<int> kept_curves, std::ofstream& f, bool log_progress = false);

	void per_curve_lambda();
};
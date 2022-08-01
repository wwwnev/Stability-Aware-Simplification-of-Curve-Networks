#include "cn_greedy_solver.h"
#include "cn_eigen_solver.h"
#include <string>
#include <execution>
#include <chrono>
#include "cn_newton_solver.h"

cn_greedy_solver::cn_greedy_solver(const curve_network& cn)
{
	m_cn = cn;
	m_cn.compute_per_curve_length();
}

void cn_greedy_solver::solve_with_relax(Eigen::VectorXd& alphas, const std::string& outname, double budget, std::vector<int> keep, std::vector<int> kick)
{
	// budget
	m_cn.compute_per_curve_length();
	if (budget == 0)
		budget = m_cn.m_curve_len.dot(alphas) / 10;

	// Every iteration :
	// Update current lambda and alphas (alpha_i = 0 -> equality constraint)
	// Stop if budget is respected
	// new cn_eigen_solver object
	// Call brute force solve 1 iter method from that object (must also return current best lambda and alphas that respect the budget, if any)

	int nb_curves = m_cn.m_curves.size();
	alphas = Eigen::VectorXd::Ones(nb_curves);
	const auto start_time = std::chrono::high_resolution_clock::now();
	double cost;
	double cur_lambda = 0;
	int removed_curve = -1;
	std::vector<int> removed_curve_idcs = kick, kept_curve_idcs = keep;

	curve_network tmp_cn = m_cn;
	double qqt_coeff = 10;
	int iteration = 1;
	 
	for (int i : removed_curve_idcs) {
		tmp_cn.remove_curve(i);
	}

	std::ofstream f;
	if (outname != "") {
		f = std::ofstream("output/greedy/" + outname + "_greedy.txt");
		f << "budget : " << budget << std::endl;
	}
	while (true) {
		std::cout << "Greedy iteration : " << iteration << std::endl;
		
		// relaxing the network
		{
			
			cn_newton_solver n_solver(tmp_cn);
			printf("relaxing from greedy solver\n");
			n_solver.solve();
			printf("done: relaxing from greedy solver\n");
			n_solver.m_CN_current.update();

			
			tmp_cn = n_solver.m_CN_current;
			tmp_cn.display("after iter " + std::to_string(iteration - 1), false, false);
			tmp_cn.write_curves("output/greedy/" + outname + "_" + std::to_string(iteration - 1) + ".txt");
		}

		tmp_cn.compute_per_curve_length();
		cost = tmp_cn.m_curve_len.sum();

		//std::cout << "current cost : " << cost << " (budget : " << budget << ")" << std::endl;
		if (cost <= budget) {
			break;
		}
		cn_eigen_solver tmp_es(tmp_cn, qqt_coeff);

		if (outname != "") {
			f << "ITERATION " << iteration << std::endl;
			f << "qqt weight : " << tmp_es.m_qqt_weight << " (from smallest eigenvalue of complete set : " << tmp_es.m_qqt_weight / qqt_coeff << ")" << std::endl;
			f << "current cost : " << cost << std::endl;
		}

		tmp_es.greedy_solve_one_iter(budget, cur_lambda, removed_curve, kept_curve_idcs, f, true); // I want lambda, the curve to remove, the best lambda that respects the budget and the removed curve.

		std::cout << "\tremoving curve " << removed_curve << std::endl;

		removed_curve_idcs.push_back(tmp_cn.m_curves[removed_curve]);

		if (outname != "") {
			f << "----" << std::endl;
			f << "this iteration, removing curve : " << tmp_cn.m_curves[removed_curve] << std::endl;
			f << "with lambda " << cur_lambda << std::endl;
			Eigen::VectorXi tmp = Eigen::VectorXi::Ones(nb_curves);
			for (int rc : removed_curve_idcs) {
				tmp(rc) = 0;
			}

			std::string s = "a = {";
			for (int i = 0; i < nb_curves; i++) {
				s.append(std::to_string(tmp(i)));
				if (i != nb_curves - 1)
					s.append(", ");
			}
			s.append("};");
			f << s << std::endl;
			const auto mid_time = std::chrono::high_resolution_clock::now();
			f << "time : " << std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(mid_time - start_time).count() << std::endl;
			f << "===========" << std::endl << std::endl;
		}
		tmp_cn.remove_curve(tmp_cn.m_curves[removed_curve]); // example : we are in 5th iteration, the first 4 removed curves 0, 1, 2, 3. Now we remove curve 4, this is gonna be m_curves[0]. since the one step greedy solve uses tmp_cn, it will not know that it is 4, but we do.

		iteration++;
	}


	for (int rc : removed_curve_idcs) {
		alphas(rc) = 0;
	}

	//std::cout << "best lambda " << cur_lambda << std::endl;
	//std::cout << "with alphas " << alphas.transpose() << std::endl;
	const auto end_time = std::chrono::high_resolution_clock::now();
	if (outname != "") {
		f << "final" << std::endl;
		std::string s = "a = {";
		for (int i = 0; i < nb_curves; i++) {
			s.append(std::to_string((int)alphas(i)));
			if (i != nb_curves - 1)
				s.append(", ");
		}
		s.append("};");
		f << s << std::endl;
		f << "lambda: " << cur_lambda << std::endl;

		f << "time : " << std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end_time - start_time).count() << std::endl;
		f.close();
	}
	m_lambda = cur_lambda;
	return;
}
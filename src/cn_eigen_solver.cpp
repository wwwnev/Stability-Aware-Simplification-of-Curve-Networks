#include "cn_eigen_solver.h"

#include <chrono>
#include <future>
#include <thread>
#include <stdlib.h>
#include "cn_newton_solver.h"

cn_eigen_solver::cn_eigen_solver(curve_network& cn, double qqt_weight_coeff)
{
	//std::cout << "cn_eigen_solver constructor" << std::endl;
	m_cn = cn;
	m_cn.compute_per_curve_hessian_yes_no();
	m_cn.compute_per_curve_length();
	m_n = m_cn.m_V.rows() * 3;
	m_nb_curves = m_cn.m_curves.size();

	m_nb_curves = m_cn.m_curves.size();
	m_qqt_weight_coefficient = qqt_weight_coeff;
	compute_qqt();
	compute_qqt_intersections();
	compute_identity_regularizers();
	compute_betas_map();

	Eigen::VectorXd alphas = Eigen::VectorXd::Ones(m_nb_curves);
	SpMat H(m_n, m_n);
	
	layer_hessians_fast(alphas, H);
	H.prune(1.0, 1e-9);
	
	std::cout << "\tcomputing weight" << std::endl;
	m_qqt_weight = compute_smallest_eigen_value(H,4) * qqt_weight_coeff;
	std::cout << "\tqqt weight : " << m_qqt_weight << " (from eigen value "<<m_qqt_weight / qqt_weight_coeff <<")" << std::endl;
	
	m_budget = m_cn.m_curve_len.dot(alphas);

}

// we dont use this anymore, see compute_identity_regularizers() for replacement
void cn_eigen_solver::compute_qqt()
{
	//std::cout << "Computing qqt matrices" << std::endl;

	for (int c_idx = 0; c_idx < m_cn.m_curves.size(); c_idx++) {

		curve c = m_cn.m_curve_set.m_all_curves[m_cn.m_curves[c_idx]];
		c.compute_segments();
		int n = c.m_V.rows();
		int n3 = n * 3;

		Eigen::SparseMatrix<double, Eigen::ColMajor> kernel(m_cn.m_V.rows() * 3, 0);

		for (int s_idx = 0; s_idx < c.m_segments.size(); s_idx++) {

			std::vector<int> s = c.m_segments[s_idx];
			Eigen::SparseMatrix<double, Eigen::ColMajor> k(m_cn.m_V.rows() * 3, 0);
			Eigen::MatrixXd rot_center(1, 3);

			if (s_idx == 0) {
				k.conservativeResize(m_cn.m_V.rows() * 3, 3);
				rot_center = c.m_V.row(s[0]);
			}
			else if (s_idx == c.m_segments.size() - 1) {
				k.conservativeResize(m_cn.m_V.rows() * 3, 3);
				rot_center = c.m_V.row(s[s.size() - 1]);
			}
			else {
				rot_center.setZero();
				k.conservativeResize(m_cn.m_V.rows() * 3, 6);
				// Translations
				for (auto i : s) {
					auto it = m_cn.m_V_correspondance[c_idx].find(i);
					int big_idx = it->second;
					k.coeffRef(big_idx * 3, 0) = 1;
					k.coeffRef(big_idx * 3 + 1, 1) = 1;
					k.coeffRef(big_idx * 3 + 2, 2) = 1;
				}
			}
			
			// Rotations
			int ncol = k.cols();

			//rot_center.setZero();
			for (auto i : s) { 

				auto it = m_cn.m_V_correspondance[c_idx].find(i);

				int big_idx = it->second;

				// xRot
				k.coeffRef(big_idx * 3 + 1, ncol - 3) = - (c.m_V(i, 2) - rot_center(0, 2));
				k.coeffRef(big_idx * 3 + 2, ncol - 3) = (c.m_V(i, 1) - rot_center(0, 1));
				// yRot
				k.coeffRef(big_idx * 3, ncol - 2) = (c.m_V(i, 2) - rot_center(0, 2));
				k.coeffRef(big_idx * 3 + 2, ncol - 2) = - (c.m_V(i, 0) - rot_center(0, 0));
				// zRot
				k.coeffRef(big_idx * 3, ncol - 1) = - (c.m_V(i, 1) - rot_center(0, 1));
				k.coeffRef(big_idx * 3 + 1, ncol - 1) = (c.m_V(i, 0) - rot_center(0, 0));
			}
			
			k.conservativeResize(k.rows(), rank_mat(k));

			SpMat kernel_tmp(kernel.rows(), kernel.cols() + k.cols());

			for (int col = 0; col < kernel.cols(); col++) {
				kernel_tmp.col(col) = kernel.col(col);
			}

			for (int col = 0; col < k.cols(); col++) {
				kernel_tmp.col(kernel.cols() + col) = k.col(col);
			}
			kernel = kernel_tmp;
		}

		SpMat q;
		Q_basis_mat(kernel, q);
		m_qqt.push_back(q * q.transpose());
	}
		
	//std::cout << "Computing qqt matrices : done" << std::endl;

	return;
}

void cn_eigen_solver::compute_qqt_intersections()
{
	//std::cout << "Computing qqt intersection matrices. Nb of intersections : " << m_cn.m_intersection_set.m_all_intersections.size()<< std::endl;
	m_cn.update_IP();
	int inter_idx = 0;
	for (auto inter : m_cn.m_intersection_set.m_all_intersections) {

		for (int i = 0; i < m_cn.m_V.rows(); i++) {

			if ((m_cn.m_V.row(i) - inter.m_pos.transpose()).norm() < 1e-8) {

				Eigen::SparseMatrix<double, Eigen::ColMajor> k(m_cn.m_V.rows() * 3, 3); // kernel for that intersection of those 2+ curves,

				k.coeffRef(i * 3, 0) = 1;
				k.coeffRef(i * 3 + 1, 1) = 1;
				k.coeffRef(i * 3 + 2, 2) = 1;

				m_qqt_intersection.push_back(k * k.transpose());
			}
		}
		inter_idx++;
	}
	//std::cout << "size of m qqt intersection : " << m_qqt_intersection.size() << std::endl;
	//std::cout << "Computing qqt intersection matrices : done" << std::endl;

	return;
}

// this is obsolete, we dont use betas_map anymore
void cn_eigen_solver::compute_betas_map()
{
	//std::cout << "Computing betas map" << std::endl;

	int betas_count = 0;

	for (int i = 0; i < m_cn.m_intersection_set.m_all_intersections.size(); i++) {
		auto a = m_cn.m_intersection_set.m_all_intersections[i].m_curve_ids;
		int tmp_a, tmp_b, count = 0;
		for (auto c_idx : a) {
			if (count == 0) {
				tmp_a = c_idx;
				count++;
			}
			else if (count == 1) {
				tmp_b = c_idx;
				count++;
			}
		}
		std::pair<int,int> tmp_pair = std::make_pair(tmp_a, tmp_b);
		m_betas_map[i] = tmp_pair;
	}

	return;

	Eigen::VectorXd layer_indices = Eigen::VectorXd::Zero(m_cn.m_V.rows()); // if layer_indices[i] == 2 -> there's intersection at vertex i


	// First check find 2 curves intersecting
	for (int c1_idx = 0; c1_idx < m_cn.m_curves.size(); c1_idx++) {

		layer_indices.setZero();
		for (const auto& it : m_cn.m_V_correspondance[c1_idx]) {
			layer_indices(it.second) += 1;
		}
		for (int c2_idx = 0; c2_idx < m_cn.m_curves.size(); c2_idx++) {
			if (c2_idx <= c1_idx) {
				// c1_idx is always smaller than c2_idx in m_betas_map, defining betas for c2_idx >= c1_idx would be a waste
				continue;
			}

			Eigen::VectorXd layer_indices_copy = layer_indices;

			for (const auto& it : m_cn.m_V_correspondance[c2_idx]) {
				layer_indices_copy(it.second) += 1;
			}

			bool no_intersection = true;
			for (int i = 0; i < layer_indices_copy.size(); i++) {
				if (layer_indices_copy(i) > 1) {
					// as soon as we detect 1 intersection between 2 curves, a beta is needed, no matter the nb of intersections.
					auto tmp_pair = std::make_pair(c1_idx, c2_idx);
					m_betas_map[betas_count] = tmp_pair;
					betas_count++;
					break;
				}
			}
		}
	}

	//std::cout << "Computing betas map : done" << std::endl;

	return;
}

void cn_eigen_solver::compute_identity_regularizers()
{
	//std::cout << std::endl<< "Computing identity regularizer matrices" << std::endl;

	for (int c_idx = 0; c_idx < m_cn.m_curves.size(); c_idx++) {
		curve c = m_cn.m_curve_set.m_all_curves[m_cn.m_curves[c_idx]];
		c.compute_segments();
		int n = c.m_V.rows();
		int n3 = n * 3;
		int sum = 0;
		Eigen::SparseMatrix<double, Eigen::ColMajor> tmp_id_reg(m_cn.m_V.rows() * 3, m_cn.m_V.rows() * 3);
		tmp_id_reg.setZero();
		for (int s_idx = 0; s_idx < c.m_segments.size(); s_idx++) {
			std::vector<int> s = c.m_segments[s_idx];
			for (auto vertex_idx : s) {
				int cn_corresponding_v_idx = m_cn.m_V_correspondance[c_idx].find(vertex_idx)->second;
				for (int i = 0; i < 3; i++) {
					sum++;
					tmp_id_reg.coeffRef(cn_corresponding_v_idx * 3 + i, cn_corresponding_v_idx * 3 + i) = 1;
				}
			}
		}
		m_identity_reg.push_back(tmp_id_reg);
	}

	//std::cout << "Computing identity regularizer matrices : done" << std::endl;

	return;
}

void cn_eigen_solver::layer_hessians(const Eigen::VectorXd& alphas, SpMat& H, double qqt_w, double w)
{
	H.conservativeResize(m_n, m_n);
	H.setZero();
	
	for (int c1_idx = 0; c1_idx < m_cn.m_curves.size(); c1_idx++) {

		H += w * (m_cn.m_H_no_bend[c1_idx] + m_cn.m_H_no_stretch[c1_idx])
			+ alphas(c1_idx) * (m_cn.m_H_yes_bend[c1_idx] + m_cn.m_H_yes_stretch[c1_idx])
			+ m_cn.m_H_penalty[c1_idx]
			+ m_qqt_weight * (1 - alphas(c1_idx)) * m_identity_reg[c1_idx];
	}

	int beta;

	for (int i = 0; i < m_cn.m_intersection_set.m_all_intersections.size(); i++) {
		bool intersection_alive = false;

		for (auto it = m_cn.m_intersection_set.m_all_intersections[i].m_curve_ids.begin(); it != m_cn.m_intersection_set.m_all_intersections[i].m_curve_ids.end(); it++) {
			int c_idx = -1;

			for (int j = 0; j < m_cn.m_curves.size(); j++) {
				if (m_cn.m_curves[j] == *it) {
					c_idx = j;
					break;
				}
			}

			if (c_idx == -1) {
				std::cerr << "BAD : We are manipulating a curve that is not currently in the curve network (previously removed)" << std::endl;
			}

			if (alphas(c_idx) == 1) {
				intersection_alive = true;
				break;
			}
		}

		intersection_alive == true ? beta = 0 : beta = 1;
		H += beta * m_qqt_weight * m_qqt_intersection[i];
	}
	H.prune(1.0, 1e-9);
	//std::cout << "Layer hessians : done" << std::endl;
	return;
}

void cn_eigen_solver::layer_hessians_fast(const Eigen::VectorXd& alphas, SpMat& H, double qqt_w, double w)
{
	H.conservativeResize(m_n, m_n);
	H.setZero();

	for (int c1_idx = 0; c1_idx < m_cn.m_curves.size(); c1_idx++) {
		H += w * (m_cn.m_H_no_bend[c1_idx] + m_cn.m_H_no_stretch[c1_idx])
			+ alphas(c1_idx) * (m_cn.m_H_yes_bend[c1_idx] + m_cn.m_H_yes_stretch[c1_idx])
			+ m_cn.m_H_penalty[c1_idx]
			+ m_qqt_weight * (1 - alphas(c1_idx)) * m_identity_reg[c1_idx];
	}
	H.prune(1.0, 1e-9);
	//std::cout << "Layer hessians : done" << std::endl;
	return;
}

std::vector<double> cn_eigen_solver::brute_force_solve_v2(double greedy_lambda, double budget)
{
	int nb_curves = m_cn.m_curves.size();
	SpMat H(m_n, m_n);

	double lambda;
	Eigen::VectorXd all_lambdas = Eigen::VectorXd::Ones(int(pow(2, nb_curves))) * -1;

	double max_lambda = 0;
	Eigen::VectorXd max_alphas;

	const auto start_time = std::chrono::high_resolution_clock::now();
	std::vector<unsigned long long int> v(pow(2, nb_curves));
	std::iota(std::begin(v), std::end(v), 0); // fills v with int such that v[i] = i

	bool slow = false;

	if (slow) {
		for (int i = 0; i < pow(2, nb_curves); i++) {
			Eigen::VectorXd a;
			SpMat H(m_n, m_n);
			convert_int_to_alphas(i, nb_curves, a);
			if (a.sum() < 1 || m_cn.m_curve_len.dot(a) > budget) {
				continue;
			}
			layer_hessians(a, H, m_qqt_weight);
			all_lambdas(i) = compute_smallest_eigen_value(H);
			std::cout << a.transpose() << " -> " << all_lambdas(i) << " (cost : " << m_cn.m_curve_len.dot(a) << ")" << std::endl;
		}
	}
	else {
		std::for_each(std::execution::par_unseq, std::begin(v), std::end(v), [&budget, &all_lambdas, &nb_curves, &start_time, this](int i) { // im afraid capturing "this" might be costly ...
			Eigen::VectorXd a;
			SpMat H(m_n, m_n);
			convert_int_to_alphas(i, nb_curves, a);
			if (a.sum() < 1 || m_cn.m_curve_len.dot(a) > budget) {
				return;
			}
			layer_hessians(a, H, m_qqt_weight);
			all_lambdas(i) = compute_smallest_eigen_value(H);
			return;
			});
		
	}

	max_lambda = all_lambdas.maxCoeff();
	int count = 0;
	double avg_lambda = 0;
	std::vector<double> lambdas_vec;
	for (int i = 0; i < all_lambdas.size(); i++) {
		if (all_lambdas(i) > 0) {
			count++;
			avg_lambda += all_lambdas(i);
			lambdas_vec.push_back(all_lambdas(i));
		}
	}
	avg_lambda /= count;

	std::sort(lambdas_vec.begin(), lambdas_vec.end());
	count = 0;

	double percentile = 0.0;
	if (abs(greedy_lambda - max_lambda) < 1e-8) {
		percentile = 1.0;
	}
	else {
		for (int i = 0; i < lambdas_vec.size(); i++) {

			if (lambdas_vec[i] > greedy_lambda) {
				break;
			}
			count++;
		}
		percentile = (double)count / (double)lambdas_vec.size();
	}

	std::vector<double> results;
	results.push_back(max_lambda);
	results.push_back(avg_lambda);
	results.push_back(percentile);
	return results;
}

void cn_eigen_solver::generate_legal_combination(double budget, const std::vector<int>& kick, const std::vector<int>& keep, Eigen::VectorXd& c) {
	int nb_curves = c.size();
	
	while (true) {

		double tmp_nb_curves = rand() % nb_curves;
		double p = (double)tmp_nb_curves / (double)nb_curves;
		double tmp2;
		for (int j = 0; j < nb_curves; j++) {
			tmp2 = (double)rand() / RAND_MAX;
			if (tmp2 > p)
				c(j) = 1;
			else
				c(j) = 0;
		}
		for (int j : keep) { // just go over the random values we just set
			c(j) = 1;
		}
		for (int j : kick) {
			c(j) = 0;
		}
		if (c.sum() < 2) {
			//std::cout << "too many zeroes !" << std::endl;
			continue;
		}
		if (m_cn.m_curve_len.dot(c) >= budget) {
			//std::cout << "too big" << std::endl;
			continue;
		}
		break;
	}
}

void cn_eigen_solver::parallel_for_brute_force(std::future<void> futureObj, Eigen::VectorXd& all_lambdas, const Eigen::MatrixXd& all_combs, int greedy_budget, const std::vector<int>& greedy_kick, const std::vector<int>& keep) {
	int nb_combinations = all_lambdas.size();
	int nb_curves = m_cn.m_curves.size();

	std::vector<unsigned long long int> v(nb_combinations);
	std::iota(std::begin(v), std::end(v), 0);

	int cur_seed = 0;

	std::for_each(std::execution::par_unseq, std::begin(v), std::end(v), [&futureObj, &all_combs, &all_lambdas, &keep, &greedy_kick, &greedy_budget, &nb_curves, this](int i) {

		if (futureObj.wait_for(std::chrono::milliseconds(1)) != std::future_status::timeout) {
			return;
		}

		SpMat H(m_n, m_n);
		layer_hessians(all_combs.row(i), H, m_qqt_weight);
		double lambda = compute_smallest_eigen_value(H);
		all_lambdas(i) = lambda;

		return;
		});
}

std::vector<double> cn_eigen_solver::random_sample_brute_force_v2(double greedy_lambda, double greedy_time, double greedy_budget, Eigen::VectorXd& max_alphas, std::vector<int> greedy_kick, const std::string& outname, std::vector<int> keep)
{

	double maximum, average, median;
	m_cn.compute_per_curve_length();
	int nb_curves = m_cn.m_curves.size();
	SpMat H(m_n, m_n);
	std::vector<Eigen::VectorXd> done_combs;
	Eigen::VectorXd alphas(nb_curves), lambdas(0), all_lambdas;

	int count = 0, tries = 0, tmp, tmp_nb_curves;
	long long int i, max_i = pow(2, nb_curves);
	double lambda;

	std::vector<double> lambdas_smth(2000, 0.0);

	std::ofstream f;
	if (outname != "") {
		f = std::ofstream("output/greedy/" + outname + "_random.txt");
	}
	double coeff_budget = greedy_budget / m_cn.m_curve_len.sum();
	double p = coeff_budget * ((double)(nb_curves - keep.size()) / (double)nb_curves);
	double tmp2;

	const auto start_time = std::chrono::high_resolution_clock::now();
	const auto cur_time = std::chrono::high_resolution_clock::now();
	double now_time = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(cur_time - start_time).count();

	double time_for_one_iter = 0;
	Eigen::MatrixXd all_combs;

	{
		auto one_iter_start = std::chrono::high_resolution_clock::now();

		std::cout << "keep : ";
		for (auto e : keep) std::cout << e << ", ";
		std::cout << std::endl;

		generate_legal_combination(greedy_budget, greedy_kick, keep, alphas);

		layer_hessians(alphas, H, m_qqt_weight);
		lambda = compute_smallest_eigen_value(H);

		time_for_one_iter = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(std::chrono::high_resolution_clock::now() - one_iter_start).count();

		int nb_of_combinations = (int)ceil(greedy_time / time_for_one_iter);

		nb_of_combinations *= 100; // let's make it big enough
		all_lambdas = Eigen::VectorXd::Zero(nb_of_combinations);
		all_combs = Eigen::MatrixXd::Zero(nb_of_combinations, nb_curves);

		for (int i = 0; i < all_combs.rows(); i++) {
			Eigen::VectorXd a = Eigen::VectorXd::Zero(nb_curves);
			generate_legal_combination(greedy_budget, greedy_kick, keep, a);
			all_combs.row(i) = a;
		}

		{
			// Create a std::promise object
			std::promise<void> exitSignal;

			//Fetch std::future object associated with promise
			std::future<void> futureObj = exitSignal.get_future();

			// Starting Thread & move the future object in lambda function by reference
			auto time_a = std::chrono::high_resolution_clock::now();

			std::thread th([&] {this->parallel_for_brute_force(std::move(futureObj), all_lambdas, all_combs, greedy_budget, greedy_kick, keep); });

			//Wait for 10 sec
			int time_left = (int)(greedy_time - std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(std::chrono::high_resolution_clock::now() - start_time).count());
			time_left = time_left < greedy_time / 20 ? (int)ceil(greedy_time / 20) : time_left;
			std::cout << "(in random) greedy time : " << greedy_time << std::endl;
			std::cout << "(in random) time left for random " << time_left << std::endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(time_left));

			//std::cout << "Asking Thread to Stop" << std::endl;

			//Set the value in promise
			exitSignal.set_value();

			//Wait for thread to join
			th.join();

			auto time_b = std::chrono::high_resolution_clock::now();
			double time_used = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(time_b - time_a).count();
			
			std::cout << "(inside random) time after par unseq : " << time_used << std::endl;

			double eliminate_ratio = time_left / time_used;

			if (eliminate_ratio < 0.99) { // ran for too long -> we computed too many combinations, let's rebuild list
				int nnz = 0;
				Eigen::VectorXd new_all_lambdas = Eigen::VectorXd::Zero(all_lambdas.size());
				Eigen::MatrixXd new_all_combs = Eigen::MatrixXd::Zero(all_combs.rows(), all_combs.cols());
				for (int i = 0; i < all_lambdas.size(); i++) {
					if (all_lambdas(i) > 0) {
						//std::cout << all_lambdas(i) << std::endl;
						new_all_lambdas(nnz) = all_lambdas(i);
						new_all_combs.row(nnz) = all_combs.row(i);
						nnz++;
					}
				}

				int to_keep = (int)ceil(nnz * eliminate_ratio), l_count = 0;

				all_lambdas = new_all_lambdas.head(to_keep);
				all_combs = new_all_combs.topRows(to_keep);
			}
		}
	}
	//std::cout << "done computing the random lambdas" << std::endl;

	if (all_lambdas.size() == 0 || all_lambdas.sum() < 1e-9) {
		return { 0,0,0,0 };
	}

	std::vector<double> results;

	lambdas = all_lambdas;

	double max_lambda = 0.0;
	int max_lambda_idx = -1;

	for (int i = 0; i < lambdas.size(); i++) {
		if (lambdas(i) > max_lambda) {
			max_lambda = lambdas(i);
			max_lambda_idx = i;
		}
	}

	max_alphas = all_combs.row(max_lambda_idx);
	std::cout << "MAX ALPHAS : " << max_alphas.transpose();

	curve_network tmp_cn = m_cn;
	for (int i = 0; i < max_alphas.size(); i++) {
		if (max_alphas(i) < 0.1)
			tmp_cn.remove_curve(i);
	}
	m_cn.display("max lambda before relax ", false);

	cn_newton_solver n_solver(tmp_cn);
	n_solver.solve();
	cn_eigen_solver e_solver(n_solver.m_CN_current, 10);
	std::cout << "max lambda before relax " << max_lambda << std::endl;
	max_lambda = e_solver.m_qqt_weight / 10;
	std::cout << "max lambda after relax " << max_lambda << std::endl;
	
	n_solver.m_CN_current.display("max lambda after relax ",false);
	n_solver.m_CN_current.write_curves("output/greedy/relaxed/soft_stadium_random_relaxed.txt");

	results.push_back(max_lambda);
	results.push_back(lambdas.mean());

	std::vector<double> lambdas_vec;
	for (int i = 0; i < lambdas.size(); i++) {
		lambdas_vec.push_back(lambdas(i));
	}

	std::sort(lambdas_vec.begin(), lambdas_vec.end());
	results.push_back(lambdas_vec[(int)lambdas_vec.size() / 2]);
	int count2 = 0;
	double percentile = 0.0;

	if (abs(greedy_lambda - results[0]) < 1e-8) {
		percentile = 1.0;
	}
	else {
		for (int i = 0; i < lambdas_vec.size(); i++) {
			//std::cout << lambdas_vec[i] << std::endl;

			if (lambdas_vec[i] > greedy_lambda) {
				break;
			}
			count2++;
		}
		percentile = (double)count2 / (double)lambdas_vec.size();
	}
	results.push_back(percentile);

	std::cout << "greedy lambda : " << greedy_lambda << std::endl;
	std::cout << "max " << results[0] << std::endl;
	std::cout << "avg " << results[1] << std::endl;
	std::cout << "med " << results[2] << std::endl;
	std::cout << "pct " << results[3] << std::endl;

	if (outname != "") {
		f << "greedy lambda : " << greedy_lambda << std::endl;
		f << "max " << results[0] << std::endl;
		f << "avg " << results[1] << std::endl;
		f << "med " << results[2] << std::endl;
		f << "pct " << results[3] << std::endl;
		f << "a = {";
		for (int i = 0; i < max_alphas.size() - 1; i++) {
			f << max_alphas(i) << ", ";
		}
		f <<max_alphas(max_alphas.size() - 1)<< "};" << std::endl;

		f.close();
	}
	return results;
}

void cn_eigen_solver::greedy_solve_one_iter(double budget, double& cur_lambda, int& removed_curve, std::vector<int> kept_curves, std::ofstream& f, bool log_progress) {

	int nb_curves = m_cn.m_curves.size();
	SpMat H(m_n, m_n);
	Eigen::VectorXd cur_alphas = Eigen::VectorXd::Ones(m_nb_curves);

	SpMat tmp_H(m_n, m_n);

	layer_hessians_fast(cur_alphas, tmp_H, m_qqt_weight);
	double tmp_lambda_full_alphas = compute_smallest_eigen_value(tmp_H);

	const auto start_time = std::chrono::high_resolution_clock::now();
	std::vector<unsigned long long int> v(m_nb_curves);
	std::iota(std::begin(v), std::end(v), 0); // fills v with int such that v[i] = i

	int iter = 1;

	double best_backup_lambda = 0.0;
	Eigen::VectorXd best_backup_alphas;

	Eigen::VectorXd iter_lambdas = Eigen::VectorXd::Zero(m_nb_curves);

	bool slow = false;

	if (slow) {
		for (int i = 0; i < m_nb_curves; i++) {
			//if (cur_alphas(i) != 0) {
			if (std::find(kept_curves.begin(), kept_curves.end(), m_cn.m_curves[i]) == kept_curves.end()) {
				Eigen::VectorXd a = cur_alphas;
				a(i) = 0;
				SpMat H(m_n, m_n);
				layer_hessians_fast(a, H, m_qqt_weight);
				double tmp_lambda = compute_smallest_eigen_value(H,1);
				std::cout << "remove i = " << i << "/" << m_nb_curves << " (curve " << m_cn.m_curves[i] << ") ... lambda = " << tmp_lambda << std::endl;
				iter_lambdas(i) = tmp_lambda;
				
			}
		}
	}
	else {
		std::for_each(std::execution::par_unseq, std::begin(v), std::end(v), [&kept_curves, &budget, &cur_alphas, &nb_curves, &start_time, &iter_lambdas, this](int i) {

			//if (cur_alphas(i) != 0) {
			if (std::find(kept_curves.begin(), kept_curves.end(), m_cn.m_curves[i]) == kept_curves.end()) {
				Eigen::VectorXd a = cur_alphas;
				a(i) = 0;
				SpMat H(m_n, m_n);
				layer_hessians_fast(a, H, m_qqt_weight);
				double tmp_lambda = compute_smallest_eigen_value(H);
				iter_lambdas(i) = tmp_lambda;			
			}
			return;
		});
	}
	Eigen::VectorXd a = cur_alphas;

	cur_lambda = 0;
	for (int i = 0; i < m_nb_curves; i++) {
		if (iter_lambdas(i) > cur_lambda) {
			cur_lambda = iter_lambdas(i);
			removed_curve = i;
		}
	}

	a(removed_curve) = 0;

	const auto end_time = std::chrono::high_resolution_clock::now();
		
	double cost = m_cn.m_curve_len.dot(a);

	if (log_progress) {
		double set_lambda = m_qqt_weight / m_qqt_weight_coefficient;
		bool monotonic = true;
		for (int i = 0; i < m_cn.m_curves.size(); i++) {
			std::string s, s2;
			s = iter_lambdas(i) > set_lambda ? "Not monotonic, ratio : " + std::to_string(iter_lambdas(i) / set_lambda) : "";
			s2 = i == removed_curve ? s2 = "** " : (iter_lambdas(i) == cur_lambda ? s2 = "*  " : "   ");
			f << s2 <<"- c" << m_cn.m_curves[i] << " -> " << std::setprecision(10)<< iter_lambdas(i) <<". "<< s<< std::endl;
		}
	}

	//std::cout << "best lambda found : " << cur_lambda << std::endl << "with " << a.transpose() << std::endl;

	return;

}

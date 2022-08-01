#include "cn_newton_solver.h"
#include "aux_tools.h"

cn_newton_solver::cn_newton_solver(curve_network cn) {
	m_CN_initial = cn;
	m_CN_current = cn;
}

void cn_newton_solver::init() {

}

double cn_newton_solver::compute_energy_and_gradient(const Eigen::VectorXd& x, Eigen::VectorXd& grad, bool compute_gradient)
{
	grad = Eigen::VectorXd::Zero(x.size());

	double fx = 0.0;

	for (int c_idx = 0; c_idx < m_CN_current.m_curves.size(); c_idx++) {
		const curve& c = m_CN_current.m_curve_set.m_all_curves[m_CN_current.m_curves[c_idx]];
		const auto& mvc = m_CN_current.m_V_correspondance[c_idx];
		if (gc_do_stretch) {
			for (int i = 0; i < c.m_E.rows(); i++) {
				std::array<int, 2> ps_idx = { mvc.at(c.m_E(i, 0)), mvc.at(c.m_E(i, 1)) };
				Eigen::Vector3d p0, p1;
				p0 << x[3 * ps_idx[0]], x[3 * ps_idx[0] + 1], x[3 * ps_idx[0] + 2];
				p1 << x[3 * ps_idx[1]], x[3 * ps_idx[1] + 1], x[3 * ps_idx[1] + 2];
				double e;				
				Codegen::computeEdgeStretchingEnergy(c.m_stiffness_modifier * gc_stretching_stiffness, p0, p1, c.m_init_e_len[i], e);
				fx += e;

				if (compute_gradient) {

					Eigen::Matrix<double, 6, 1> mini_grad_stretch0;
					Codegen::computeEdgeStretchingEnergyGradient(c.m_stiffness_modifier * gc_stretching_stiffness, p0, p1, c.m_init_e_len[i], mini_grad_stretch0);

					for (int j = 0; j < 6; j++) {
						int j_flr = floor(j / 3);
						int j_rmdr = j % 3;
						grad(3 * ps_idx[j_flr] + j_rmdr) += mini_grad_stretch0(j);
					}
				}
			}
		}

		if (gc_do_bend) {
			int nb_bends = c.m_closed ? c.m_E.rows() : c.m_E.rows() - 1;
			for (int i = 0; i < nb_bends; i++) {
				std::array<int, 3> ps_idx;
				double n0, n1;

				if (i == c.m_E.rows() - 1) {
					ps_idx = { mvc.at(c.m_E(i, 0)), mvc.at(c.m_E(i, 1)), mvc.at(c.m_E(0, 1)) };
					n0 = c.m_init_e_len[i];
					n1 = c.m_init_e_len[0];
				}
				else {
					ps_idx = { mvc.at(c.m_E(i, 0)), mvc.at(c.m_E(i, 1)), mvc.at(c.m_E(i + 1, 1)) };
					n0 = c.m_init_e_len[i];
					n1 = c.m_init_e_len[i+1];
				}
				Eigen::Vector3d p0, p1, p2;
				p0 << x[3 * ps_idx[0]], x[3 * ps_idx[0] + 1], x[3 * ps_idx[0] + 2];
				p1 << x[3 * ps_idx[1]], x[3 * ps_idx[1] + 1], x[3 * ps_idx[1] + 2];
				p2 << x[3 * ps_idx[2]], x[3 * ps_idx[2] + 1], x[3 * ps_idx[2] + 2];


				double e;
				Codegen::computeStraightBendingEnergy(c.m_stiffness_modifier * gc_bending_stiffness, p0, p1, p2, n0, n1, e);
				fx += e;
				//fxs[c_idx] += e;

				if (compute_gradient) {

					Eigen::Matrix<double, 9, 1> mini_grad;
					Codegen::computeStraightBendingEnergyGradient(c.m_stiffness_modifier * gc_bending_stiffness, p0, p1, p2, n0, n1, mini_grad);

					// fit stuff into the gradient
					for (int j = 0; j < 9; j++) {
						int j_flr = floor(j / 3);
						int j_rmdr = j % 3;
						grad(3 * ps_idx[j_flr] + j_rmdr) += mini_grad(j) /*+ mini_grad_stretch[j]*/;

					}
				}
			}
		}
		if (gc_fix_endpoints && !c.m_closed) {
			int idx;
			double dist;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 2; j++) {

					idx = mvc.find(c.m_E(0, j))->second;

					dist = x(3 * idx + i) - c.m_V_initial(c.m_E(0, j), i);
					fx += gc_endpoint_penalty_mu * dist * dist; // E = mu * ||x - c||^2 -> dE/dx = 2 * mu (x - c)
					
					if (compute_gradient)
						grad(3 * idx + i) += 2 * gc_endpoint_penalty_mu * dist;

					idx = mvc.find(c.m_E(c.m_E.rows() - 1, j))->second;

					dist = x(3 * idx + i) - c.m_V_initial(c.m_E(c.m_E.rows() - 1, j), i);
					fx += gc_endpoint_penalty_mu * dist * dist;

					if (compute_gradient)
						grad(3 * idx + i) += 2 * gc_endpoint_penalty_mu * dist;
				}

			}

		}

	}
	return fx;
}

double cn_newton_solver::compute_energy(const Eigen::VectorXd& x){
	Eigen::VectorXd tmp;
	return compute_energy_and_gradient(x, tmp, false);
}

void cn_newton_solver::solve(const std::string& s, int max_i, double epsilon)
{
	std::cout << "Starting newton method descent" << std::endl;
	double fx, step, p, c, fx_w_step, suff_decrease;
	if (epsilon < 0) // will be -1 if not provided
		epsilon = 1e-2;
	if (max_i == 0)
		max_i = 1000;
	int cur_i = 0;

	Eigen::VectorXd grad, tmp_grad, x;
	double last_step_size = 1;
	double last_fx = 0;
	int use_reg_count = 0;

	int n = m_CN_initial.m_V.rows() * 3;
	SpMat reg = SpMat(n, n);
	double cee = 1e-6;

	bool do_grad_fd = false;
	bool do_hessian_fd = false;
	bool do_grad_descent = false;

	std::vector<Eigen::Triplet<double>> triplets;
	for (int i = 0; i < n; i++)
		triplets.push_back({ i,i,1 });
	reg.setFromTriplets(triplets.begin(), triplets.end());

	do {
		mat_to_vector(m_CN_current.m_V, x);

		// Computing the descent direction
		m_CN_current.generate_elastic_energy_hessian(false);
		grad.setZero();
		fx = compute_energy_and_gradient(x, grad);
		std::cout << "\titer "<<cur_i<<"\tfx : " << fx <<"\tgrad norm "<<grad.norm()<< std::endl;

		// every iteration m_V saved, will try using them to reconstruct surface.
		if (s != "") {
			std::string ss = s + std::to_string(cur_i) + ".txt";
			m_CN_current.write_m_V(ss);
		}

		if (grad.norm() < epsilon) {
			auto tmp_grad = grad;
			auto grad_abs = tmp_grad.cwiseAbs();
			int l0_norm = 0;
			for (int i = 0; i < grad.size(); i++) {
				if (fabs(grad(i)) > 1e-9)
					l0_norm++;
			}

			std::cout << "\tGradient small : "<<grad.norm()<<". Stopping after " << cur_i << " iterations" << std::endl;

			if (s != "") {
				std::string ss = s + std::to_string(cur_i) + ".txt";
			}
			break;
		}

		SpMat H = m_CN_current.m_H;
		
		int k = 0;
		Eigen::VectorXd delta_x;
		if (do_grad_descent)
			delta_x = grad;
		else {
			auto init_H = H;
			do {
				Eigen::SimplicialLDLT<SpMat> chol(H);

				if (chol.info() != Eigen::Success) {
					// decomposition failed
					std::cerr << "newton method : decomposition failed at iter: " << cur_i << std::endl;
					return;
				}

				delta_x = chol.solve(grad);

				if ((H*delta_x - grad).norm() > 0.01)
					std::cout << "RESIDUAL BIG : " << (H * delta_x - grad).norm() << std::endl;

				k += 1;
				H = m_CN_current.m_H + cee * pow(2, k) * reg;

			} while (-delta_x.dot(grad) > -1e-10);
			if (k > 1) {
				std::cout << "\t\t k = " << k << std::endl;
				use_reg_count++;
			}
		}

		
		bool using_grad = false;
		if (-delta_x.transpose().dot(grad) > 0) {
			//std::cout << "descent direction dot gradient is positive -> using -grad as descent direction for this iteration instead of newton" << std::endl;
			delta_x = grad;
			use_reg_count++;
			using_grad = true;
		}

		// Choosing step size using sufficient decrease and backtracking
		bool do_naive_step_selection = true;
		if (do_naive_step_selection) {
			step = 1;
			p = 0.90; // how do i make it vary. contraction factor for sufficient decrease and backtracking
			c = 1e-4;

			std::vector<double> fx_w_steps;
			std::vector<double> suff_decreases;

			fx_w_step = compute_energy(x -step * delta_x);
			suff_decrease = fx + c * step * -grad.transpose().dot(delta_x);

			int tries = 0;
			while (fx_w_step > suff_decrease) {
				step = step * p;
				fx_w_step = compute_energy(x - step * delta_x);

				suff_decrease = fx + c * step * -grad.transpose().dot(delta_x);

				tries++;
			}
			last_step_size = step;

		}

		// Updating m_CN_current
		std::cout << "\tstep : " << step << std::endl;
		Eigen::MatrixXd move;
		vector_to_mat(step * delta_x, delta_x.size()/3, 3, move);
		m_CN_current.m_V -= move;
		m_CN_current.update_curves();

		if (cur_i % 500 == 1 && false) {
			if (s != "") {
				std::string ss = s + "_m_V_" + std::to_string(cur_i) + ".txt";
			}
		}

		cur_i += 1;

	} 	while (cur_i < max_i);

	m_CN_current.update_IP();
	m_CN_current.update();

	if (s != "") {
		std::string ss = s + "_m_V_" + std::to_string(cur_i) + ".txt";
	}

	//std::cout << "Ending newton method descent" << std::endl;

	return;
}

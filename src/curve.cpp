#include "curve.h"

curve::curve()
{
	m_E = Eigen::MatrixXi(0, 2);
	m_V = Eigen::MatrixXd(0, 3);
	m_stiffness = gc_bending_stiffness;
	m_stiffness_modifier = 1.0;
}

// New edges are appended to the bottom of m_E
bool curve::split_edge(int e, std::vector<Eigen::Vector3d> inter_pts) {

	if (inter_pts.size() == 0) return false;

	auto compare_dist = [] (Eigen::Vector3d a, Eigen::Vector3d b) {
		return (a.norm() < b.norm());
	};

	// we want to sort inter_pts by distance to the first vertex of edge e
	Eigen::Vector3d p0 = m_V.row(m_E(e, 0)), p1 = m_V.row(m_E(e, 1));

	for (auto it = inter_pts.begin(); it != inter_pts.end(); ) {
		if (((*it - p0).norm() < 1e-8) || ((*it - p1).norm() < 1e-8)) {
			//std::cout << "erasing a splitting point because it's on the vertices of the edge" << std::endl;
			it = inter_pts.erase(it);
		}
		else {
			++it;
		}
	}

	if (inter_pts.size() == 0)
		return true;

	for (int i = 0; i < inter_pts.size(); i++) {
		inter_pts[i] -= p0;
	}

	std::sort(inter_pts.begin(), inter_pts.end(), compare_dist);

	for (int i = 0; i < inter_pts.size(); i++) {
		inter_pts[i] += p0;
	}

	if ((inter_pts[0].transpose() - m_V.row(m_E(0, 0))).norm() < 1e-3) {
		m_V.row(m_E(0, 0)) = inter_pts[0].transpose();
		inter_pts.erase(inter_pts.begin());
	}
	else if ((inter_pts[0].transpose() - m_V.row(m_E(m_E.rows() - 1, 1))).norm() < 1e-3) {
		m_V.row(m_E(m_E.rows() - 1, 1)) = inter_pts[0].transpose();
		inter_pts.erase(inter_pts.begin());
	}
	if (inter_pts.size() == 0) {
		return true;
	}

	if ((inter_pts[inter_pts.size() - 1].transpose() - m_V.row(m_E(0, 0))).norm() < 1e-3) {
		m_V.row(m_E(0, 0)) = inter_pts[inter_pts.size() - 1].transpose();
		inter_pts.erase(inter_pts.end() - 1);
	}
	else if ((inter_pts[inter_pts.size() - 1].transpose() - m_V.row(m_E(m_E.rows() - 1, 1))).norm() < 1e-3) {
		m_V.row(m_E(m_E.rows() - 1, 1)) = inter_pts[inter_pts.size() - 1].transpose();
		inter_pts.erase(inter_pts.end() - 1);
	}
	if (inter_pts.size() == 0) {
		return true;
	}

	// replacing the edge by the first segment
	int p_idx = add_vertex(inter_pts[0]);
	m_E(e, 1) = p_idx;
	m_inter_pts.insert(p_idx);

	// only goes into this loop when there is 2 or more splitting points
	for (int i = 0; i < inter_pts.size() - 1; i++) {
		std::pair<int, int> a = add_edge(inter_pts[i], inter_pts[i + 1], e + i);
		m_inter_pts.insert(std::get<1>(a));
	}

	// adding one last edge
	add_edge(inter_pts[inter_pts.size() - 1], p1, e + inter_pts.size() - 1);

	//std::cout << "Splitting edges of a curve : done" << std::endl;

	return true;
}

void curve::compute_segments()
{
	// want to build the m_segments member. std::vector<Eigen::MatrixXi>. The eigen mat are pieces of m_E.

	m_segments.clear();
	
	std::vector<int> one_seg;
	bool first_vertex_is_ip = false;


	for (auto ip : m_inter_pts) {
		if (m_E(0, 0) == ip) {
			first_vertex_is_ip = true;
			break;
		}
	}

	for (int row_idx = 0; row_idx < m_E.rows(); row_idx++) {
		bool found_ip = false;
		for (auto ip : m_inter_pts) {
			if (m_E(row_idx, 0) == ip) {
				found_ip = true;
				if (one_seg.size()>0) {
					m_segments.push_back(one_seg);
					one_seg.clear();
				}
				break;
			}
		}
		if (!found_ip) {
			one_seg.push_back(m_E(row_idx, 0));
		}
	}
	bool last_vertex_is_ip = false;
	for (auto ip : m_inter_pts) {
		if (m_E(m_E.rows() - 1, 1) == ip) {

			last_vertex_is_ip = true;
			break;
		}
	}
	if (!last_vertex_is_ip)
		one_seg.push_back(m_E(m_E.rows() - 1, 1));

	if (one_seg.size() > 0) // in case the last vertex of the curve is an intersection vertex, we dont want to push an empty segment.
		m_segments.push_back(one_seg);

	
	if (m_closed && !first_vertex_is_ip) { // we want to join last and first segment. Currently, this never happens because first and last vertex have an intersection in cn.m_intersections (from cn::add_curve)
		std::vector<int> a = m_segments[m_segments.size() - 1];

		for (int i = 0; i < m_segments[0].size(); i++) {
			if (m_segments[0][i] != a[a.size() - 1])
				a.push_back(m_segments[0][i]);
		}

		m_segments[0] = a;

		m_segments.pop_back();
	}
	return;

}

std::pair<int,int> curve::add_edge(Eigen::MatrixXd v0, Eigen::MatrixXd v1)
{
	int idx_v0 = add_vertex(v0);
	int idx_v1 = add_vertex(v1);
	m_E.conservativeResize(m_E.rows() + 1, 2);
	m_E(m_E.rows() - 1, 0) = idx_v0;
	m_E(m_E.rows() - 1, 1) = idx_v1;
	return std::make_pair(idx_v0, idx_v1);
}

// takes v0 and v1, the 2 end vertices of the edge, and the edge index of the precedent edge -> tells us where to insert the edge in m_E.
std::pair<int, int> curve::add_edge(Eigen::MatrixXd v0, Eigen::MatrixXd v1, int e_precedent)
{
	int idx_v0 = add_vertex(v0);
	int idx_v1 = add_vertex(v1);

	//inserting edge in the "middle" of m_E
	Eigen::MatrixXi newE(m_E.rows() + 1, 2);

	newE.topRows(e_precedent + 1) = m_E.topRows(e_precedent + 1);
	newE.bottomRows(m_E.rows() - e_precedent - 1) = m_E.bottomRows(m_E.rows() - e_precedent - 1);
	newE(e_precedent + 1, 0) = idx_v0;
	newE(e_precedent + 1, 1) = idx_v1;
	m_E = newE;

	return std::make_pair(idx_v0, idx_v1);
}

int curve::add_vertex(Eigen::MatrixXd v)
{
	double epsilon = 1e-8;
	for (int i = 0; i < m_V.rows(); i++) {
		if ((v.transpose() - m_V.row(i)).norm() >= 0 - epsilon && (v.transpose() - m_V.row(i)).norm() <= 0 + epsilon) {
			return i;
		}
	}

	m_V.conservativeResize(m_V.rows() + 1, 3);
	int idx = m_V.rows() - 1;
	
	if (v.rows() == 3 && v.cols() == 1)
		m_V.row(idx) = v.transpose();
	else // we assume this function is called only with 3x1 or 1x3 matrices / vectors
		m_V.row(idx) = v;

	return idx;
}

bool curve::cut(Eigen::Vector3d p)
{
	int e_split_idx, v_split_idx;
	double len1 = 0, len2 = 0;
	bool do_one = true;
	for (int e_idx = 0; e_idx < m_E.rows(); e_idx++) {
		if (do_one)
			len1 += (m_V.row(m_E(e_idx, 0)) - m_V.row(m_E(e_idx, 1))).norm();
		else
			len2 += (m_V.row(m_E(e_idx, 0)) - m_V.row(m_E(e_idx, 1))).norm();
		if ((p.transpose() - m_V.row(m_E(e_idx, 1))).norm() < 1e-8) {
			v_split_idx = m_E(e_idx, 1);
			e_split_idx = e_idx;
			do_one = false;
		}
		

	}
	std::cout << " len 1 " << len1 << " len 2 " << len2 << std::endl;
	Eigen::MatrixXd new_V(0, 3);
	Eigen::MatrixXi new_E(0, 2);

	if (len1 > len2) {
		new_V.conservativeResize(new_V.rows() + 1, 3);
		new_V.row(new_V.rows() - 1) = m_V.row(m_E(0, 0));

		for (int e_idx = 0; e_idx <= e_split_idx; e_idx++) {
			new_V.conservativeResize(new_V.rows() + 1, 3);
			new_E.conservativeResize(new_E.rows() + 1, 2);
			new_V.row(new_V.rows() - 1) = m_V.row(m_E(e_idx, 1));
			new_E(new_E.rows() - 1, 0) = new_V.rows() - 2;
			new_E(new_E.rows() - 1, 1) = new_V.rows() - 1;
		}
	}
	else {
		new_V.conservativeResize(new_V.rows() + 1, 3);

		new_V.row(new_V.rows() - 1) = p;
		for (int e_idx = e_split_idx + 1; e_idx < m_E.rows(); e_idx++) {
			new_V.conservativeResize(new_V.rows() + 1, 3);
			new_E.conservativeResize(new_E.rows() + 1, 2);
			new_V.row(new_V.rows() - 1) = m_V.row(m_E(e_idx, 1));
			new_E(new_E.rows() - 1, 0) = new_V.rows() - 2;
			new_E(new_E.rows() - 1, 1) = new_V.rows() - 1;
			
		}
	}
	m_V = new_V;
	m_E = new_E;
	return true;
}

double curve::get_length() {
	double l = 0.0;
	for (int i = 0; i < m_E.rows(); i++) {
		l += (m_V.row(m_E(i, 1)) - m_V.row(m_E(i, 0))).norm();
	}
	return l;
}

// edges in m_E must be ordered (e.g. 0 -- 1, 1 -- 10, 10 -- 3, ...) for disambiguation of rod direction at high valence vertices
// Only works on open curves. It should be easy to modify to take into consideration if the curve is closed or not.
// But it is a bit weird to think that a closed curve was once a straight rod. (different topology)
void curve::compute_full_elastic_energy_hessian() {
	int n = m_V.rows() * 3;
	m_H = SpMat(n,n);
	m_H.setZero();
	std::vector<Eigen::Triplet<double>> triplets;
	if (gc_do_bend) {

		int nb_bends = m_closed? m_E.rows() : m_E.rows() - 1;

		for (int e0_idx = 0; e0_idx < nb_bends; e0_idx++) {
			std::vector<int> ps_idx;

			double n0, n1;
			if (e0_idx == m_E.rows() - 1) {
				ps_idx.push_back(m_E(e0_idx, 0));
				ps_idx.push_back(m_E(e0_idx, 1));
				ps_idx.push_back(m_E(0, 1));
				n0 = m_init_e_len[e0_idx];
				n1 = m_init_e_len[0];
			}
			else {
				ps_idx.push_back(m_E(e0_idx, 0));
				ps_idx.push_back(m_E(e0_idx, 1));
				ps_idx.push_back(m_E(e0_idx + 1, 1));
				n0 = m_init_e_len[e0_idx];
				n1 = m_init_e_len[e0_idx + 1];
			}

			Eigen::Vector3d p0 = m_V.row(ps_idx[0]), p1 = m_V.row(ps_idx[1]), p2 = m_V.row(ps_idx[2]);

			Eigen::Matrix<double, 9, 9> h;
			
			Codegen::computeStraightBendingEnergyHessian(m_stiffness_modifier * gc_bending_stiffness, p0, p1, p2, n0, n1, h);
			
			// now let's fit this 9x9 hessian into the hessian of the energy of the whole curve

			for (int i = 0; i < h.rows(); i++) {
				for (int j = 0; j < h.cols(); j++) {
					int i_flr = floor(i / 3), j_flr = floor(j / 3);
					int i_rmdr = i % 3, j_rmdr = j % 3;

					// indices in the big hessian. row and col used. Symmetry of the matrix is handled by the 2 for loops
					int row_ind = 3 * ps_idx[i_flr] + i_rmdr;
					int col_ind = 3 * ps_idx[j_flr] + j_rmdr;

					triplets.push_back({ row_ind, col_ind, h(i,j) });
				}
			}
			
		}
	}
	
	if (gc_do_stretch) {
		for (int e0_idx = 0; e0_idx < m_E.rows(); e0_idx++) {

			std::vector<int> ps_idx;
			ps_idx.push_back(m_E(e0_idx, 0));
			ps_idx.push_back(m_E(e0_idx, 1));

			Eigen::Vector3d P0 = m_V_initial.row(ps_idx[0]), P1 = m_V_initial.row(ps_idx[1]);
			Eigen::Vector3d p0 = m_V.row(ps_idx[0]), p1 = m_V.row(ps_idx[1]);

			Eigen::Matrix<double, 6, 6> mini_hessian_stretch0;
			
			Codegen::computeEdgeStretchingEnergyHessian(m_stiffness_modifier * gc_stretching_stiffness, p0, p1, m_init_e_len[e0_idx], mini_hessian_stretch0);

			// now let's fit this 6x6 hessian into the hessian of the energy of the whole curve
			for (int i = 0; i < mini_hessian_stretch0.rows(); i++) {
				for (int j = 0; j < mini_hessian_stretch0.cols(); j++) {
					int i_flr = floor(i / 3), j_flr = floor(j / 3);
					int i_rmdr = i % 3, j_rmdr = j % 3;

					int row_ind = 3 * ps_idx[i_flr] + i_rmdr;
					int col_ind = 3 * ps_idx[j_flr] + j_rmdr;

					triplets.push_back({ row_ind, col_ind, mini_hessian_stretch0(i,j) });
				}
			}
		}
	}
	
	double mu = gc_endpoint_penalty_mu;
	if (gc_fix_endpoints && !m_closed) {
		for (int i = 0; i < 3; i++) {
			int row_idx = 3 * m_E(0, 0) + i;
			triplets.push_back({ row_idx, row_idx, 2 * mu});

			row_idx = 3 * m_E(0, 1) + i;
			triplets.push_back({ row_idx, row_idx, 2 * mu }); // enforcing tangents

			row_idx = 3 * m_E(m_E.rows() - 1, 0) + i;
			triplets.push_back({ row_idx, row_idx, 2 * mu }); // enforcing tangents

			row_idx = 3 * m_E(m_E.rows() - 1, 1) + i;
			triplets.push_back({ row_idx, row_idx, 2 * mu });
		}
		
	}
	m_H.setFromTriplets(triplets.begin(), triplets.end());
	
	return;
}

void curve::compute_yes_no_bending_hessian(SpMat& yes, SpMat& no, const std::set<int>& ip)
{
	int n = m_V.rows() * 3;
	yes = SpMat(n, n);
	yes.setZero();
	no = SpMat(n, n);
	no.setZero();

	std::vector<Eigen::Triplet<double>> triplets_yes, triplets_no;

	int nb_bends = m_closed ? m_E.rows() : m_E.rows() - 1;
	for (int e0_idx = 0; e0_idx <nb_bends; e0_idx++) {
		
		std::vector<int> ps_idx;
		double n0, n1;
		if (e0_idx == m_E.rows() - 1) {
			ps_idx.push_back(m_E(e0_idx, 0));
			ps_idx.push_back(m_E(e0_idx, 1));
			ps_idx.push_back(m_E(0, 1));
			n0 = m_init_e_len[e0_idx];
			n1 = m_init_e_len[0];
		}
		else {
			ps_idx.push_back(m_E(e0_idx, 0));
			ps_idx.push_back(m_E(e0_idx, 1));
			ps_idx.push_back(m_E(e0_idx + 1, 1));
			n0 = m_init_e_len[e0_idx];
			n1 = m_init_e_len[e0_idx + 1];
		}

		Eigen::Vector3d p0 = m_V.row(ps_idx[0]), p1 = m_V.row(ps_idx[1]), p2 = m_V.row(ps_idx[2]);

		Eigen::Matrix<double, 9, 9> h;

		Codegen::computeStraightBendingEnergyHessian(m_stiffness_modifier * gc_bending_stiffness, p0, p1, p2, n0, n1, h);

		// now let's fit this 9x9 hessian into the hessian of the energy of the whole curve
		bool yes_alpha = (ip.find(ps_idx[0]) != ip.end() || ip.find(ps_idx[1]) != ip.end() || ip.find(ps_idx[2]) != ip.end());

		for (int i = 0; i < h.rows(); i++) {
			for (int j = 0; j < h.cols(); j++) {
				int i_flr = floor(i / 3), j_flr = floor(j / 3);
				int i_rmdr = i % 3, j_rmdr = j % 3;

				// indices in the big hessian. row and col used. Symmetry of the matrix is handled by the 2 for loops
				int row_ind = 3 * ps_idx[i_flr] + i_rmdr;
				int col_ind = 3 * ps_idx[j_flr] + j_rmdr;

				if (yes_alpha) {
					triplets_yes.push_back({row_ind, col_ind, h(i,j)});
				}
				else {
					triplets_no.push_back({ row_ind, col_ind, h(i,j) });
				}
			}
		}
	}
	yes.setFromTriplets(triplets_yes.begin(), triplets_yes.end());
	no.setFromTriplets(triplets_no.begin(), triplets_no.end());

	return;
}

void curve::compute_yes_no_stretching_hessian(SpMat& yes, SpMat& no, const std::set<int>& ip)
{
	int n = m_V.rows() * 3;
	yes = SpMat(n, n);
	yes.setZero();
	no = SpMat(n, n);
	no.setZero();

	std::vector<Eigen::Triplet<double>> triplets_yes, triplets_no;

	for (int e0_idx = 0; e0_idx < m_E.rows(); e0_idx++) {

		std::vector<int> ps_idx;
		ps_idx.push_back(m_E(e0_idx, 0));
		ps_idx.push_back(m_E(e0_idx, 1));

		Eigen::Vector3d p0 = m_V.row(ps_idx[0]), p1 = m_V.row(ps_idx[1]);

		Eigen::Matrix<double, 6, 6> mini_hessian_stretch0;
		Codegen::computeEdgeStretchingEnergyHessian(m_stiffness_modifier * gc_stretching_stiffness, p0, p1, m_init_e_len[e0_idx], mini_hessian_stretch0);

		// now let's fit this 6x6 hessian into the hessian of the energy of the whole curve
		bool yes_alpha = (ip.find(ps_idx[0]) != ip.end() || ip.find(ps_idx[1]) != ip.end());

		for (int i = 0; i < mini_hessian_stretch0.rows(); i++) {
			for (int j = 0; j < mini_hessian_stretch0.cols(); j++) {
				int i_flr = floor(i / 3), j_flr = floor(j / 3);
				int i_rmdr = i % 3, j_rmdr = j % 3;

				int row_ind = 3 * ps_idx[i_flr] + i_rmdr;
				int col_ind = 3 * ps_idx[j_flr] + j_rmdr;

				if (yes_alpha) {
					triplets_yes.push_back({ row_ind, col_ind, mini_hessian_stretch0(i,j) });
				}
				else {
					triplets_no.push_back({ row_ind, col_ind, mini_hessian_stretch0(i,j) });
				}
			}
		}
	}
	yes.setFromTriplets(triplets_yes.begin(), triplets_yes.end());
	no.setFromTriplets(triplets_no.begin(), triplets_no.end());
	return;
}
////////
double curve::compute_init_lengths()
{
	m_init_e_len.resize(m_E.rows(), 0);
	double l = 0.0;
	double el = 0.0;
	for (int i = 0; i < m_E.rows(); i++) {
		el = (m_V.row(m_E(i, 1)) - m_V.row(m_E(i, 0))).norm();
		l += el;
		m_init_e_len[i] = el;
	}
	return l;
}

void curve::compute_penalty_hessian(SpMat& penalty)
{
	int n = m_V.rows() * 3;
	penalty = SpMat(n, n);
	penalty.setZero();

	std::vector<Eigen::Triplet<double>> triplets;
	double mu = gc_endpoint_penalty_mu;
	int row_idx;
	if (gc_fix_endpoints && !m_closed) {
		for (int i = 0; i < 3; i++) {
			row_idx = 3 * m_E(0, 0) + i;
			triplets.push_back({ row_idx, row_idx, 2 * mu });
			row_idx = 3 * m_E(0, 1) + i;
			//triplets.push_back({ row_idx, row_idx, 2 * mu }); // enforcing tangents

			row_idx = 3 * m_E(m_E.rows() - 1, 0) + i;
			//triplets.push_back({ row_idx, row_idx, 2 * mu }); // enforcing tangents
			row_idx = 3 * m_E(m_E.rows() - 1, 1) + i;
			triplets.push_back({ row_idx, row_idx, 2 * mu });
		}
	}
	penalty.setFromTriplets(triplets.begin(), triplets.end());

	return;
}

void curve::set_m_V_initial() {
	m_V_initial = m_V;
}

void curve::resample(double l) {
	//std::cout << "Resampling curve" << m_V.rows() << " " << m_E.rows()<< std::endl;

	

	//std::cout << "ips : ";
	//for (auto e : m_inter_pts) std::cout << e << ", ";
	//std::cout << std::endl;

	// note : segment = piece of a curve, delimited by an end point + intersection point, or 2 intersection pts. for this fct, we are talking specifically about the smallest possible segments.

	// copying the old curve
	std::set<int> m_inter_pts_copy = m_inter_pts;
	Eigen::MatrixXd m_V_copy = m_V;
	Eigen::MatrixXi m_E_copy = m_E;

	// clearing the old curve
	m_inter_pts.clear();
	m_V.resize(0, 3);
	m_E.resize(0, 2);

	int v_start = m_E_copy(0, 0);
	int e_start = 0;

	int v_end;
	int e_end;

	// the length of every edges
	std::vector<double> edge_lengths;
	for (int e_idx = 0; e_idx < m_E_copy.rows(); e_idx++) {
		double len = (m_V_copy.row(m_E_copy(e_idx, 0)) - m_V_copy.row(m_E_copy(e_idx, 1))).norm();
		edge_lengths.push_back(len);
	}

	Eigen::VectorXd prev_pt, next_pt;
	prev_pt = m_V_copy.row(v_start);

	int seg = 0;
	bool finished = false;

	while (!finished) {

		// Lets find the end of this curve segment (intersection point or end of curve)
		bool found = false;
		for (int e_idx = e_start; e_idx < m_E_copy.rows(); e_idx++) {
			int next_v_idx = m_E_copy(e_idx, 1);
			for (int inter_pt_idx : m_inter_pts_copy) {
				if (next_v_idx == inter_pt_idx) {
					// reached an intersection point
					v_end = next_v_idx;
					e_end = e_idx;
					found = true;
					break;
				}
			}
			if (found) break;
		}
		// In the case where there were no intersection points on the remaining segment, we reached the end of the curve
		if (!found) {
			e_end = m_E_copy.rows() - 1;
			v_end = m_E_copy(e_end, 1);
			finished = true;
		}

		// now we know all the edges that must be split. 

		// What's the current length of the segment
		double total_length = 0.0;
		for (int e_idx = e_start; e_idx <= e_end; e_idx++) {
			total_length += edge_lengths[e_idx];
		}

		// how many edges do we want in the segment ?
		int num_new_edges = round(total_length / l);
		num_new_edges = num_new_edges < 6 ? 6 : num_new_edges;
		double new_l = total_length / num_new_edges;
		double epsilon = 1e-8;
		std::tuple<int, int> ae;
		for (int i = 0; i < num_new_edges; i++) {
			double curr_length = 0.0;
			double needed_length = (i + 1) * new_l;

			// looking for the point at needed_length distance from v_start
			int j;
			for (j = e_start; j <= e_end; j++) {
				if (curr_length + edge_lengths[j] + epsilon > needed_length) {
					//std::cout << "hellooooo" << std::endl;
					break; // we know that the new point is gonna be on edge j
				}
				else curr_length += edge_lengths[j];
			}
			double alpha = (needed_length - curr_length) / edge_lengths[j];
			next_pt = m_V_copy.row(m_E_copy(j, 0)) * (1 - alpha) + m_V_copy.row(m_E_copy(j, 1)) * (alpha);
			ae = add_edge(prev_pt, next_pt);
			prev_pt = next_pt;
		}

		next_pt = m_V_copy.row(v_end);

		// band aid fix, might never actually pass this condition ...
		if ((prev_pt - next_pt).norm() > epsilon)
			ae = add_edge(prev_pt, next_pt);

		prev_pt = next_pt;

		if (!finished) {
			m_inter_pts.insert(std::get<1>(ae));
			e_start = e_end + 1;
			v_start = m_E_copy(e_start, 0);
		}
		
	}
	//set_m_V_initial();
	//std::cout << "Resampling curve : done" <<m_V.rows()<<" "<<m_E.rows()<< std::endl;
}

void curve::remove_0_edges() {

	std::vector<int> v_to_remove;

	for (int i = 0; i < m_E.rows(); i++) {
		double e_len = (m_V.row(m_E(i, 1)) - m_V.row(m_E(i, 0))).norm();
		if (e_len < 1e-8) {
			if (m_E(i, 0) == m_E(i, 1))
				continue;
			// removing vertex with highest idx
			int to_remove, to_keep;
			to_remove = m_E(i, 1) > m_E(i, 0) ? m_E(i, 1) : m_E(i, 0);
			to_keep = m_E(i, 1) > m_E(i, 0) ? m_E(i, 0) : m_E(i, 1);

			Eigen::MatrixXd new_m_V(m_V.rows() - 1, 3);
			int count = 0;
			for (int j = 0; j < m_V.rows(); j++) {
				if (j == to_remove) {
					continue;
				}
				new_m_V.row(count) = m_V.row(j);
				count++;
			}
			m_V = new_m_V;
			for (int j = 0; j < m_E.rows(); j++) {
				for (int k = 0; k < 2; k++) {
					if (m_E(j, k) == to_remove) {
						m_E(j, k) = to_keep;
					}
					// need to shift the index of everything in m_E by -1, if the index was bigger than to_remove
					if (m_E(j, k) > to_remove) {
						m_E(j, k) -= 1;
					}
				}
			}
			
		}
	}
	Eigen::MatrixXi new_m_E(0, 2);
	for (int i = 0; i < m_E.rows(); i++) {
		if (m_E(i, 0) == m_E(i, 1)) {
			continue;
		}
		new_m_E.conservativeResize(new_m_E.rows() + 1, 2);
		new_m_E.row(new_m_E.rows() - 1) = m_E.row(i);
	}
	m_E = new_m_E;
	return;
}

void curve::resample(double l, const Eigen::MatrixXd& intersection_vertices)
{
	//std::cout << "Resampling curve" << m_V.rows() << " " << m_E.rows() << std::endl;
	remove_0_edges();
	bool any_nan = false;

	for (int i = 0; i < m_V.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			if (isnan(m_V(i, j))) {
				any_nan = true;
				break;
			}
		}
		if (any_nan)
			break;
	}
	if (any_nan) {
		std::cout << "found some NaNs" << std::endl;
	}

	//std::cout << "ips : ";
	//std::cout << intersection_vertices << std::endl;
	//std::cout << "m_E : " <<std::endl<< m_E << std::endl;
	//std::cout << "m_V : " << std::endl<<m_V << std::endl;
	
	// note : segment = piece of a curve, delimited by an end point + intersection point, or 2 intersection pts. for this fct, we are talking specifically about the smallest possible segments.

	// copying the old curve
	// std::set<int> m_inter_pts_copy = m_inter_pts;
	Eigen::MatrixXd m_V_copy = m_V;
	Eigen::MatrixXi m_E_copy = m_E;

	// clearing the old curve
	m_inter_pts.clear();
	m_V.resize(0, 3);
	m_E.resize(0, 2);

	int v_start = m_E_copy(0, 0);
	int e_start = 0;

	int v_end;
	int e_end;

	// the length of every edges
	std::vector<double> edge_lengths;
	for (int e_idx = 0; e_idx < m_E_copy.rows(); e_idx++) {
		double len = (m_V_copy.row(m_E_copy(e_idx, 0)) - m_V_copy.row(m_E_copy(e_idx, 1))).norm();
		edge_lengths.push_back(len);
	}

	Eigen::VectorXd prev_pt, next_pt;
	prev_pt = m_V_copy.row(v_start);

	int seg = 0;
	bool finished = false;

	while (!finished) {

		// Lets find the end of this curve segment (intersection point or end of curve)
		bool found = false;
		Eigen::Vector3d v1, v2;
		for (int e_idx = e_start; e_idx < m_E_copy.rows(); e_idx++) {
			int next_v_idx = m_E_copy(e_idx, 1);
			v1 = m_V_copy.row(next_v_idx);
			for (int inter_pt_idx = 0; inter_pt_idx < intersection_vertices.rows(); inter_pt_idx++) {
				v2 = intersection_vertices.row(inter_pt_idx);
				if ((v1 - v2).norm() < 1e-9) {
					v_end = next_v_idx;
					e_end = e_idx;
					found = true;
					break;
				}
			}
			//for (int inter_pt_idx : m_inter_pts_copy) {
			//	if (next_v_idx == inter_pt_idx) {
			//		// reached an intersection point
			//		v_end = next_v_idx;
			//		e_end = e_idx;
			//		found = true;
			//		break;
			//	}
			//}
			if (found) break;
		}
		// In the case where there were no intersection points on the remaining segment, we reached the end of the curve
		if (!found) {
			e_end = m_E_copy.rows() - 1;
			v_end = m_E_copy(e_end, 1);
			finished = true;
		}
		else if (v_end == m_E_copy(m_E_copy.rows() - 1, 1)) {
			finished = true;
		}

		// now we know all the edges that must be split. 

		// What's the current length of the segment
		double total_length = 0.0;
		for (int e_idx = e_start; e_idx <= e_end; e_idx++) {
			total_length += edge_lengths[e_idx];
		}

		// how many edges do we want in the segment ?
		int num_new_edges = round(total_length / l);
		num_new_edges = num_new_edges < 2 ? 2 : num_new_edges;
		double new_l = total_length / num_new_edges;
		double epsilon = 1e-8;
		std::tuple<int, int> ae;
		for (int i = 0; i < num_new_edges; i++) {
			double curr_length = 0.0;
			double needed_length = (i + 1) * new_l;

			// looking for the point at needed_length distance from v_start
			int j;
			for (j = e_start; j <= e_end; j++) {
				if (curr_length + edge_lengths[j] + epsilon > needed_length) {
					//std::cout << "hellooooo" << std::endl;
					break; // we know that the new point is gonna be on edge j
				}
				else curr_length += edge_lengths[j];
			}
			double alpha = (needed_length - curr_length) / edge_lengths[j];
			next_pt = m_V_copy.row(m_E_copy(j, 0)) * (1 - alpha) + m_V_copy.row(m_E_copy(j, 1)) * (alpha);
			
			for (int k = 0; k < 3; k++) {
				if (isnan(next_pt(k))) {
					std::cout<<"   NAN POINT FOR RESAMPLING... new len per edge "<<new_l<<", length of segment "<<total_length<<", e_start "<<e_start<<", e_end "<<e_end<<", len e_start"<< edge_lengths[e_start]<<std::endl;
					break;
				}
			}
			ae = add_edge(prev_pt, next_pt);
			prev_pt = next_pt;
		}

		next_pt = m_V_copy.row(v_end);

		// band aid fix, might never actually pass this condition ...
		if ((prev_pt - next_pt).norm() > epsilon)
			ae = add_edge(prev_pt, next_pt);

		prev_pt = next_pt;

		if (!finished) {
			m_inter_pts.insert(std::get<1>(ae));
			e_start = e_end + 1;
			v_start = m_E_copy(e_start, 0);
		}

	}
	//set_m_V_initial();


	any_nan = false;

	for (int i = 0; i < m_V.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			if (isnan(m_V(i, j))) {
				any_nan = true;
				break;
			}
		}
		if (any_nan)
			break;
	}
	if (any_nan) {
		std::cout << "found some NaNs" << std::endl;
		std::cout << intersection_vertices << std::endl << std::endl;
	}

	//std::cout << "Resampling curve : done" << m_V.rows() << " " << m_E.rows() << std::endl;
}

void curve::replace_edge_by_vertex(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& pnew)
{
	// Do I need to worry about having multiple p1 or p2 vertices ?

	int v1_idx = -1, v2_idx = -1, tmp_idx;
	// find vertex index of p1 p2 -> v_idx1, v_idx2
	for (int i = 0; i < m_V.rows(); i++) {
		if ((m_V.row(i) - p1.transpose()).norm() < 1e-10) {
			v1_idx = i;
		}
	}
	for (int i = m_V.rows() - 1; i >= 0; i--) {
		if ((m_V.row(i) - p2.transpose()).norm() < 1e-10) {
			v2_idx = i;
		}
	}

	if (v2_idx == v1_idx) {
		std::cerr << "replace_edge_by_vertex : the vertices defining the edges provided are one and the same." << std::endl;
	}
	// we want v1_idx < v2_idx always.
	else if (!(v1_idx > v2_idx)) {
		tmp_idx = v1_idx;
		v1_idx = v2_idx;
		v2_idx = tmp_idx;
	}

	// add pnew to m_V (replace p1 by pnew in fact)
	m_V.row(v1_idx) = pnew;

	// remove p2
	Eigen::MatrixXd tmp_m_V(m_V.rows() - 1, m_V.cols());
	int count = 0;
	for (int i = 0; i < m_V.rows(); i++) {
		if (i != v2_idx) {
			tmp_m_V.row(count) = m_V.row(i);
			count++;
		}
	}
	m_V = tmp_m_V;
	
	//for (int i = 0; i < m_V.rows(); i++) {
	//	std::cout << i << " -> " << m_V.row(i) << std::endl;
	//}

	// replace every v2_idx by v1_idx
	for (int i = 0; i < m_E.rows(); i++) {
		if (m_E(i, 0) == v2_idx)
			m_E(i, 0) = v1_idx;
		else if (m_E(i, 1) == v2_idx)
			m_E(i, 1) = v1_idx;
	}
	// if an edge has m_E(i,0) == m_E(i,1), just remove it, it's useless.

	Eigen::MatrixXi tmp_m_E(0, m_E.cols());
	for (int i = 0; i < m_E.rows(); i++) {
		if (m_E(i, 0) != m_E(i, 1)) {
			//std::cout << i;
			tmp_m_E.conservativeResize(tmp_m_E.rows() + 1, tmp_m_E.cols());
			int tmp_idx = tmp_m_E.rows() - 1;
			tmp_m_E.row(tmp_idx) = m_E.row(i);
			if (tmp_m_E(tmp_idx, 0) > v2_idx)
				tmp_m_E(tmp_idx, 0) -= 1;
			if (tmp_m_E(tmp_idx, 1) > v2_idx)
				tmp_m_E(tmp_idx, 1) -= 1;
			//std::cout << i<<std::endl;
		}
	}
	m_E = tmp_m_E;
}

void curve::update_inter_pts_segments(const intersection_set& is)
{
	m_inter_pts.clear();
	
	for (intersection a : is.m_all_intersections) {
		for (int i = 0; i < m_V.rows(); i++) { 
			if ((a.m_pos - m_V.row(i).transpose()).norm() < 1e-9) {
				m_inter_pts.insert(i);
			}
		}
	}

	compute_segments();
}

void curve::close_curve()
{
	
	
	if (!check_closed()) {
		m_E.conservativeResize(m_E.rows() + 1, 2);
		m_E(m_E.rows() - 1, 0) = m_E(m_E.rows() - 2, 1);
		m_E(m_E.rows() - 1, 1) = m_E(0, 0);
		if (check_closed()) std::cout << "the curve is now closed" << std::endl;
		else std::cout << std::endl << "***** this curve cannot be closed ? investigate" << std::endl << std::endl;
	}
	else std::cout << "the curve was already closed" << std::endl;
	return;
}

bool curve::check_closed()
{
	if (m_E(0, 0) == m_E(m_E.rows() - 1, 1)) {
		m_closed = true;
	}
	else if ((m_V.row(m_E(0, 0)) - m_V.row(m_E(m_E.rows() - 1, 1))).norm() < 1e-5) {
		std::cout << std::endl << "****this is a closed curve, BUT the first and last pt of the curve do not share the same vertex idx" << std::endl<<std::endl;
		m_closed = true;
	}
	else {
		m_closed = false;
	}
	return m_closed;
}

bool curve::join(curve& c)
{
	if (c.m_V.rows() == 0)
		return false;
	if ((m_V.row(m_E(0,0)) - c.m_V.row(c.m_E(0,0))).norm() > 1e-5)
		return false;

	Eigen::MatrixXi tmp_m_E(m_E.rows() + c.m_E.rows(), 2);
	Eigen::MatrixXd tmp_m_V(m_V.rows() + c.m_V.rows() - 1, 3);

	tmp_m_V.row(0) = c.m_V.row(c.m_E(c.m_E.rows() - 1, 1));
	int v0_idx = 0, v1_idx = 1;
	int count = 0;
	for (int i = c.m_E.rows() - 1; i >= 0; i--) {
		tmp_m_V.row(count + 1) = c.m_V.row(c.m_E(i,0));
		tmp_m_E(count, 0) = v0_idx;
		tmp_m_E(count, 1) = v1_idx;
		v0_idx += 1;
		v1_idx += 1;
		count++;
	}

	for (int i = 0; i < m_E.rows(); i++) {
		tmp_m_V.row(i + c.m_V.rows()) = m_V.row(m_E(i, 1));
		tmp_m_E(count, 0) = v0_idx;
		tmp_m_E(count, 1) = v1_idx;
		v0_idx += 1;
		v1_idx += 1;
		count++;
	}
	
	m_V = tmp_m_V;
	m_E = tmp_m_E;
	return true;

}

void curve::replace_closest_vertex(Eigen::Vector3d new_pos) {
	//std::vector<double> dists;
	int min_dist_idx = -1;
	double min_dist = 0;
	for (int i = 0; i < m_V.rows(); i++) {
		double dist = (new_pos.transpose() - m_V.row(i)).norm();

		if (min_dist_idx == -1 || dist < min_dist) {
			min_dist_idx = i;
			min_dist = dist;
		}
	}

	m_V.row(min_dist_idx) = new_pos;

	return;
}

void curve::smooth_new_segments(const Eigen::Vector3d& new_vertex, const Eigen::MatrixXd& intersection_vertices)
{
	//std::cout << "Resampling curve" << m_V.rows() << " " << m_E.rows() << std::endl;
	

	//std::cout << "ips : ";
	//std::cout << intersection_vertices << std::endl;
	//std::cout << "m_E : " <<std::endl<< m_E << std::endl;
	//std::cout << "m_V : " << std::endl<<m_V << std::endl;

	// note : segment = piece of a curve, delimited by an end point + intersection point, or 2 intersection pts. for this fct, we are talking specifically about the smallest possible segments.

	// copying the old curve
	// std::set<int> m_inter_pts_copy = m_inter_pts;
	Eigen::MatrixXd m_V_copy = m_V;
	Eigen::MatrixXi m_E_copy = m_E;

	// clearing the old curve
	m_inter_pts.clear();
	m_V.resize(0, 3);
	m_E.resize(0, 2);

	int v_start = m_E_copy(0, 0);
	int e_start = 0;

	int v_end;
	int e_end;

	// the length of every edges
	std::vector<double> edge_lengths;
	for (int e_idx = 0; e_idx < m_E_copy.rows(); e_idx++) {
		double len = (m_V_copy.row(m_E_copy(e_idx, 0)) - m_V_copy.row(m_E_copy(e_idx, 1))).norm();
		edge_lengths.push_back(len);
	}

	Eigen::VectorXd prev_pt, next_pt;
	prev_pt = m_V_copy.row(v_start);

	int seg = 0;
	bool finished = false;
	std::cout << "before while " << std::endl;
	while (!finished) {
		std::cout << e_start << std::endl;
		// Lets find the end of this curve segment (intersection point or end of curve)
		bool found = false;
		Eigen::Vector3d v1, v2;
		for (int e_idx = e_start; e_idx < m_E_copy.rows(); e_idx++) {
			int next_v_idx = m_E_copy(e_idx, 1);
			v1 = m_V_copy.row(next_v_idx);
			for (int inter_pt_idx = 0; inter_pt_idx < intersection_vertices.rows(); inter_pt_idx++) {
				v2 = intersection_vertices.row(inter_pt_idx);
				if ((v1 - v2).norm() < 1e-9) {
					v_end = next_v_idx;
					e_end = e_idx;
					found = true;
					break;
				}
			}
			//for (int inter_pt_idx : m_inter_pts_copy) {
			//	if (next_v_idx == inter_pt_idx) {
			//		// reached an intersection point
			//		v_end = next_v_idx;
			//		e_end = e_idx;
			//		found = true;
			//		break;
			//	}
			//}
			if (found) break;
		}
		// In the case where there were no intersection points on the remaining segment, we reached the end of the curve
		if (!found) {
			e_end = m_E_copy.rows() - 1;
			v_end = m_E_copy(e_end, 1);
			finished = true;
		}
		else if (v_end == m_E_copy(m_E_copy.rows() - 1, 1)) {
			finished = true;
		}

		// new stuff for smoothing :
		bool start_is_new;
		if ((new_vertex.transpose() - m_V_copy.row(v_start)).norm() < 1e-8)
			start_is_new = true;
		else if ((new_vertex.transpose() - m_V_copy.row(v_end)).norm() < 1e-8)
			start_is_new = false;
		else {
			if (!finished) {
				e_start = e_end + 1;
				v_start = m_E_copy(e_start, 0);
			}
			continue;
		}
		Eigen::MatrixXd new_m_V = m_V_copy;
		Eigen::Vector3d p0, p1, p2, new_p;
		for (int e_idx = e_start; e_idx < e_end; e_idx++) {
			p0 = m_V_copy.row(m_E_copy(e_idx, 0));
			p1 = m_V_copy.row(m_E_copy(e_idx, 1)); // this vertex will be moved.
			p2 = m_V_copy.row(m_E_copy(e_idx + 1, 1));

			new_m_V.row(m_E_copy(e_idx, 1)) = (p0 + p1 + p2) / 3;
		}
		m_V_copy = new_m_V;

		next_pt = m_V_copy.row(v_end);
		prev_pt = next_pt;
		if (!finished) {
			e_start = e_end + 1;
			v_start = m_E_copy(e_start, 0);
		}

	}
	m_V = m_V_copy;
	m_E = m_E_copy;
	//set_m_V_initial();
	return;
}

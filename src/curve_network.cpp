#include <filesystem>

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include "curve_network.h"

curve_network::curve_network() {
	m_V = Eigen::MatrixXd(0, 3);
	m_E = Eigen::MatrixXi(0, 2);
	m_stiffness = gc_bending_stiffness;
}

void curve_network::add_curve(curve c) {
	//std::cout << std::endl << "ADDING CURVE " << m_curves.size() << std::endl;

	auto compare_idx = [](std::tuple<int, std::vector<Eigen::Vector3d>> a, std::tuple<int, std::vector<Eigen::Vector3d>> b) {
		return (std::get<0>(a) < std::get<0>(b));
	};

	auto remove_duplicate = [](std::vector<Eigen::Vector3d>& vec) {

		std::vector<Eigen::Vector3d> vec2;
		for (int i = 0; i < vec.size(); i++) {
			bool already_in = false;
			for (int j = 0; j < vec2.size(); j++) {
				if ((vec[i] - vec2[j]).norm() < 1e-9) {
					already_in = true;
					break;
				}
			}
			if (!already_in) {
				vec2.push_back(vec[i]);
			}
		}

		return vec2;
	};

	// test, maybe it breaks stuff ...
	c.check_closed();

	int new_curve_idx = m_curve_set.add_curve(c) - 1;
	m_curves.push_back(new_curve_idx);

	std::vector<std::tuple<int, std::vector<Eigen::Vector3d>>> pts_to_add_per_edge; // tuples (idx edge, vector of new points on the edge)
	
	double epsilon = 5e-3;
	int tmp_count = 0;

	// loop over existing curves in the CN
	for (int idx = 0; idx < m_curves.size(); idx++) {

		bool is_last;
		is_last = idx == m_curves.size() - 1 ? true : false;

		int c_idx = m_curves[idx];

		std::vector<std::tuple<int, std::vector<Eigen::Vector3d>>> pts_to_add_per_edge_old_curve;
		
		// for each segment in the new curve
		for (int e1 = 0; e1 < c.m_E.rows(); e1++) {

			// for each segment in the current curve of the curve network we want to compare to
			for (int e2 = 0; e2 < m_curve_set.m_all_curves[c_idx].m_E.rows(); e2++) {
				
				if (is_last) {
					if (e1 - 1 == e2 || e1 == e2 || e1 + 1 == e2) {
						continue;
					}
				}

				Eigen::Vector3d v0 = c.m_V.row(c.m_E(e1, 0)), v1 = c.m_V.row(c.m_E(e1, 1)), v2 = m_curve_set.m_all_curves[c_idx].m_V.row(m_curve_set.m_all_curves[c_idx].m_E(e2, 0)), v3 = m_curve_set.m_all_curves[c_idx].m_V.row(m_curve_set.m_all_curves[c_idx].m_E(e2, 1));
				bool should_intersect = false;

				if ((v1 - v0).norm() <= 1e-8 || (v3 - v2).norm() <= 1e-8) {
					continue;
				}

				if (((v1 - v0).cross(v3 - v2)).norm() <= 1e-9) {
					if (c_idx == 16 && m_curves.size() - 1 == 27) {
						tmp_count++;
					}
					continue;
				}

				if (fabs((v0 - v2).norm()) > epsilon && fabs((v0 - v3).norm()) > epsilon && fabs((v1 - v2).norm()) > epsilon && fabs((v1 - v3).norm()) > epsilon) {
					
					bool same_plane = check_coplanarity_lines(v0, v1, v2, v3);
					if (same_plane) {
						Eigen::Vector3d inter_pt;
						double t, s;
						find_line_intersection(v0, v1, v2, v3, inter_pt, t, s);
						bool intersection_exists = (s >= 0 && s <= 1 && t >= 0 && t <= 1);
						if (intersection_exists) {
							int inter_idx = m_intersection_set.add_intersection(inter_pt, c_idx, m_curves.size() - 1);
							inter_pt = m_intersection_set.m_all_intersections[inter_idx].m_pos;
							bool first_time_splitting = true;

							if (!(s > 1 - 1e-8 || s < 1e-8)) {
								for (int i = 0; i < pts_to_add_per_edge.size(); i++) {
									if (std::get<0>(pts_to_add_per_edge[i]) == e1) {
										std::get<1>(pts_to_add_per_edge[i]).push_back(inter_pt);
										first_time_splitting = false;
										break;
									}
								}
								if (first_time_splitting) {
									std::vector<Eigen::Vector3d> tmp;
									tmp.push_back(inter_pt);
									pts_to_add_per_edge.push_back(std::make_tuple(e1, tmp));
								}
							}

							// and on the old curve IF it's not the new curve (last curve of the m_curves vector
							first_time_splitting = true;
							if (!(t > 1 - 1e-8 || t < 1e-8)) {
								for (int i = 0; i < pts_to_add_per_edge_old_curve.size(); i++) {
									if (std::get<0>(pts_to_add_per_edge_old_curve[i]) == e2) {
										std::get<1>(pts_to_add_per_edge_old_curve[i]).push_back(inter_pt);
										first_time_splitting = false;
										break;
									}
								}
								if (first_time_splitting) {
									std::vector<Eigen::Vector3d> tmp;
									tmp.push_back(inter_pt);
									pts_to_add_per_edge_old_curve.push_back(std::make_tuple(e2, tmp));
								}
							}
							std::vector<int> curves_vec;
							curves_vec.push_back(m_curves.size() - 1);
							curves_vec.push_back(c_idx);
						}
					}
				}

				// they intersect at an already existing vertex.
				else if (fabs((v0 - v2).norm()) <= epsilon) {
					if (!is_last) {
						m_intersection_set.add_intersection(v2, c_idx, m_curves.size() - 1);
						m_curve_set.m_all_curves[m_curves.size() - 1].m_inter_pts.insert(c.m_E(e1, 0));
						// if the 2 points are not exactly the same, but still their dist is < epsilon, we want to replace one of them with the other, so now they are the exact same. We choose to modify the curve being added.
						m_curve_set.m_all_curves[m_curves.size() - 1].m_V.row(c.m_E(e1, 0)) = v2; // v0 = v2
					}
					m_curve_set.m_all_curves[c_idx].m_inter_pts.insert(m_curve_set.m_all_curves[c_idx].m_E(e2, 0));
				}
				else if (fabs((v1 - v2).norm()) <= epsilon) {
					if (!is_last) {
						m_intersection_set.add_intersection(v2, c_idx, m_curves.size() - 1);
						m_curve_set.m_all_curves[m_curves.size() - 1].m_inter_pts.insert(c.m_E(e1, 1));
						m_curve_set.m_all_curves[m_curves.size() - 1].m_V.row(c.m_E(e1, 1)) = v2; // v0 = v2
					}
					m_curve_set.m_all_curves[c_idx].m_inter_pts.insert(m_curve_set.m_all_curves[c_idx].m_E(e2, 0));
				}
				else if (fabs((v0 - v3).norm()) <= epsilon) {
					if (!is_last) {
						m_intersection_set.add_intersection(v3, c_idx, m_curves.size() - 1);
						m_curve_set.m_all_curves[m_curves.size() - 1].m_inter_pts.insert(c.m_E(e1, 0));
						m_curve_set.m_all_curves[m_curves.size() - 1].m_V.row(c.m_E(e1, 0)) = v3; // v0 = v3
					}
					m_curve_set.m_all_curves[c_idx].m_inter_pts.insert(m_curve_set.m_all_curves[c_idx].m_E(e2, 1));
				}
				else if (fabs((v1 - v3).norm()) <= epsilon) {
					if (!is_last) {
						m_intersection_set.add_intersection(v3, c_idx, m_curves.size() - 1);
						m_curve_set.m_all_curves[m_curves.size() - 1].m_inter_pts.insert(c.m_E(e1, 1));
						m_curve_set.m_all_curves[m_curves.size() - 1].m_V.row(c.m_E(e1, 1)) = v3; // v0 = v2
					}
					m_curve_set.m_all_curves[c_idx].m_inter_pts.insert(m_curve_set.m_all_curves[c_idx].m_E(e2, 1));
				}
			}
		}
		// would use a lambda, but m_curves is not a variable ...
		
		if (!is_last) { // if it's the very last curve, it's c itself, let's not split it here, but out of the for loop that loops over all curves in m_curve_set.m_all_curves
			std::sort(pts_to_add_per_edge_old_curve.begin(), pts_to_add_per_edge_old_curve.end(), compare_idx);

			for (int i = pts_to_add_per_edge_old_curve.size() - 1; i >= 0; i--) {
				std::tuple<int, std::vector<Eigen::Vector3d>> tupl = pts_to_add_per_edge_old_curve[i];
				int e = std::get<0>(tupl);
				std::vector<Eigen::Vector3d> l = std::get<1>(tupl);
				l = remove_duplicate(l);
				m_curve_set.m_all_curves[c_idx].split_edge(e, l);
			}
		}
	}

	std::sort(pts_to_add_per_edge.begin(), pts_to_add_per_edge.end(), compare_idx);

	for (int i = pts_to_add_per_edge.size() - 1; i >= 0; i--) {
		std::tuple<int, std::vector<Eigen::Vector3d>> tupl = pts_to_add_per_edge[i];
		int e = std::get<0>(tupl);
		std::vector<Eigen::Vector3d> l = std::get<1>(tupl);
		l = remove_duplicate(l);
		m_curve_set.m_all_curves[m_curves.size() -1 ].split_edge(e, l);
	}
	//std::cout << "Checking for intersections : done" << std::endl;

}

void curve_network::add_curve_for_display(curve c)
{
	m_curve_set.add_curve(c);
	int new_curve_idx = m_curve_set.add_curve(c) - 1;
	m_curves.push_back(new_curve_idx);
}

void curve_network::update_IP() {
	Eigen::VectorXd cur_IP = Eigen::VectorXd::Zero(m_V.rows());
	m_IP.clear();

	for (int i = 0; i < m_curves.size(); i++) {

		for (const auto& it : m_V_correspondance[i]) {
			cur_IP(it.second) += 1;
		}
	}
	for (int i = 0; i < cur_IP.size(); i++) {
		if (cur_IP(i) > 1) {
			m_IP.insert(i);
		}
	}

	// update m_intersection_set also ? PROBABLY WONT WORK FOR SELF INTERSECTION !!!
	intersection_set new_inter_set;
	for (int c1 = 0; c1 < m_curves.size() - 1; c1++) {
		int c1_idx = m_curves[c1];
		const auto& correspo1 = m_V_correspondance[c1];

		for (int c2 = c1 + 1; c2 < m_curves.size(); c2++) {
			int c2_idx = m_curves[c2];
			const auto& correspo2 = m_V_correspondance[c2];

			for (const auto& e1 : correspo1) {
				for (const auto& e2 : correspo2) {
					if (e1.second == e2.second) {
						new_inter_set.add_intersection(m_V.row(e1.second), c1_idx, c2_idx);
					}
				}
			}
		}
	}
	m_intersection_set = new_inter_set;

	int hits = 0;
	for (int i = 0; i < m_V.rows(); i++) {
		for (const auto& inter : m_intersection_set.m_all_intersections) {
			if ((m_V.row(i) - inter.m_pos.transpose()).norm() < 1e-8)
				hits++;
		}
	}
	//std::cout << "Updated IPs, checking if cn.m_V coincide with cn.m_intersection_set. Nb of hits : " << hits << std::endl;
	return;
}

std::tuple<int,int> curve_network::add_edge(Eigen::MatrixXd v0, Eigen::MatrixXd v1)
{
	int idx_v0 = add_vertex(v0);
	int idx_v1 = add_vertex(v1);
	m_E.conservativeResize(m_E.rows() + 1, 2);
	m_E(m_E.rows() - 1, 0) = idx_v0;
	m_E(m_E.rows() - 1, 1) = idx_v1;
	return std::make_tuple(idx_v0, idx_v1);
}

int curve_network::add_vertex(Eigen::MatrixXd v)
{
	double epsilon = 1e-3;
	for (int i = 0; i < m_V.rows(); i++) {
		if ((v - m_V.row(i)).norm() >= 0 - epsilon && (v - m_V.row(i)).norm() <= 0 + epsilon) {
			m_IP.insert(i); 
			return i;
		}
	}
	m_V.conservativeResize(m_V.rows() + 1, 3);
	int idx = m_V.rows() - 1;
	m_V.row(idx) = v;
	return idx;
}

void curve_network::remove_curve(int c_idx)
{
	// Before removing the curve, let's just count the number of intersection per curve, if it doesnt change by removing the curve, no need to resample.
	std::map<int, int> old_nb_intersections;
	for (auto c2_idx : m_curves) {
		std::pair<int, int> tmp_pair(c2_idx, m_intersection_set.get_nb_intersections(c2_idx));
		old_nb_intersections.insert(tmp_pair);
	}

	// We want to remove a curve from the current subset of curves (so from vector<int> m_curves).
	// Remove c_idx from m_curves


	for (int i = 0; i < m_curves.size(); i++) {
		if (m_curves[i] == c_idx) {
			m_curves.erase(m_curves.begin() + i);
			m_V_correspondance.erase(m_V_correspondance.begin() + i);
			break;
		}
	}

	// Remove c_idx from m_intersection_set
	m_intersection_set.remove_curve(c_idx);

	// Resample every curve
	for (auto c2_idx : m_curves) {
		Eigen::MatrixXd inter_vertices;
		// gather every intersection vertex position for curve with c2_idx
		int nb_inter = m_intersection_set.get_intersection_vertices(c2_idx, inter_vertices);

		if (nb_inter != old_nb_intersections.at(c2_idx)) {

			// commenting this line for now, when doing greedy algo we cannot resample because it will change the energy of the curves.

			//m_curve_set.m_all_curves[c2_idx].resample(gc_resample_length, inter_vertices);
		}
		
	}

	// Update
	//update();
	update_fast();
	set_m_V_initial();

}

// Call AFTER CURVE REMOVED FROM m_curves and m_V_correspondance. Replaces update().
void curve_network::update_fast() {

	std::map<int, int> old_to_new; // first old vertex idx, second 

	// make a new m_V matrix. more rows than it needs, we will remove the bottom rows at the end.
	Eigen::MatrixXd new_m_V = Eigen::MatrixXd::Zero(m_V.rows(), 3);
	Eigen::MatrixXi new_m_E = Eigen::MatrixXi::Zero(m_E.rows(), 2);
	std::vector<std::map<int, int>> new_m_V_correspondance;
	m_IP.clear();

	// check that m_V_correspondance has the same size as m_curves
	if (m_V_correspondance.size() != m_curves.size()) {
		std::cerr << "m_V_correspondance not same size as m_curves" << m_V_correspondance.size() << " vs " << m_curves.size() << std::endl;
		return;
	}

	// iterate over m_V_correspondance
	int v_count = 0;
	for (int c_idx = 0; c_idx < m_curves.size(); c_idx++) {
		std::map<int, int> new_correspo;
		// 
		auto correspo = m_V_correspondance[c_idx];
		std::map<int, int>::iterator it;
		for (auto v_pair : correspo) {
			it = old_to_new.find(v_pair.second);
			if (it == old_to_new.end()) { // this means we have not yet added the vertex.

				new_m_V.row(v_count) = m_V.row(v_pair.second); // adding the vertex

				new_correspo[v_pair.first] = v_count; // updating correspondance : v_idx in curve -> v_idx in curve network.

				old_to_new[v_pair.second] = v_count;

				v_count++;
			}
			else { // if the vertex is already added to new_m_V, we know it's an intersection -> add to m_IP. Probably doesnt work for self intersection. m_IP should be built from m_intersection_set
				new_correspo[v_pair.first] = it->second;
				m_IP.insert(it->second);
			}
		}
		new_m_V_correspondance.push_back(new_correspo);

	}

	m_V = new_m_V.topRows(v_count);

	int e_count = 0;
	for (int i = 0; i < m_E.rows(); i++) {

		auto it0 = old_to_new.find(m_E(i, 0));
		auto it1 = old_to_new.find(m_E(i, 1));

		if (it0 != old_to_new.end() && it1 != old_to_new.end()) { // this old edge has both of it's vertices in the new m_V, therefore we copy the edge. We can do that because we know that there will at least be one non intersection vertex between 2 intersection vertices
			new_m_E(e_count, 0) = old_to_new[m_E(i, 0)];
			new_m_E(e_count, 1) = old_to_new[m_E(i, 1)];
			e_count++;
		}
	}
	m_E = new_m_E.topRows(e_count);
	m_V_correspondance = new_m_V_correspondance;

	// We dont care about m_V_end if we dont do relaxation

	return;
}

void curve_network::update_for_display() {
	m_V.resize(0, 3);
	m_E.resize(0, 2);
	for (const curve& c : m_curve_set.m_all_curves) {
		int existing_vertices = m_V.rows();
		m_V.conservativeResize(m_V.rows() + c.m_V.rows(), 3);
		m_E.conservativeResize(m_E.rows() + c.m_E.rows(), 2);
		auto m_E_copy = c.m_E;

		m_V.bottomRows(c.m_V.rows()) = c.m_V;

		m_E_copy+= Eigen::MatrixXi::Ones(m_E_copy.rows(), 2) * existing_vertices;

		m_E.bottomRows(c.m_E.rows()) = m_E_copy;
	}
}

// Will rebuild m_V and m_E and m_V_correspondance and m_V_end. Call before doing a big operation on / with the curve network
void curve_network::update() {
	//std::cout << "UPDATING CURVE NETWORK (may take a while (too long))" << std::endl;
	m_V_correspondance.clear();
	m_V_end.clear();
	m_IP.clear();
	m_V.resize(0, 3);
	m_E.resize(0, 2);
	
	for (auto c_idx : m_curves) {
		curve c = m_curve_set.m_all_curves[c_idx];
		// append vertices
		std::map<int, int> correspondance;
		std::tuple<int, int> idcs;
		for (int i = 0; i < c.m_E.rows(); i++) {
			idcs = add_edge(c.m_V.row(c.m_E(i, 0)), c.m_V.row(c.m_E(i, 1)));//m_V_correspondance
			correspondance.insert({c.m_E(i, 0), std::get<0>(idcs)});
			if (i == 0 && !c.m_closed) {
				m_V_end.insert(std::get<0>(idcs));
				m_V_end.insert(std::get<1>(idcs)); // test : adding the before last vertices to the m_V_end member. It is used in relaxation. I hope this will help preserve the launch angle through the relaxation
			}
			else if (i == c.m_E.rows() - 1 && !c.m_closed) {
				m_V_end.insert(std::get<0>(idcs)); // test : idem
				m_V_end.insert(std::get<1>(idcs));
			}
		}
		correspondance.insert({c.m_E(c.m_E.rows() - 1, 1), std::get<1>(idcs)});
		m_V_correspondance.push_back(correspondance);
	}
	update_IP();
	//std::cout << "UPDATING CURVE NETWORK : done."<< std::endl;
	set_m_V_initial();
	return;
}

void curve_network::update_curves() {
	// with the help of m_V_correspondance, 
	// update the curve objects of m_curves with respect to m_V.
	// this function is to be called after m_V was modified.

	for (int c_idx = 0; c_idx < m_V_correspondance.size(); c_idx++) {
		for (auto cor : m_V_correspondance[c_idx]) {
			
			m_curve_set.m_all_curves[m_curves[c_idx]].m_V.row(std::get<0>(cor)) = m_V.row(std::get<1>(cor));
		}
	}
}

void curve_network::fix_m_V_end_bean_w_arcs(double z, bool use_z) {
	m_V_end.clear();

	for (int i = 0; i < m_curves.size(); i++) {
		int c_idx = m_curves[i];
		std::cout << " curve " << c_idx;
		if (m_curve_set.m_all_curves[c_idx].m_closed) {
			std::cout << " is closed" << std::endl;
		}
		else {
			std::cout << " is open" << std::endl;
			curve c = m_curve_set.m_all_curves[c_idx];

			std::map<int, int> correspondance;
			std::tuple<int, int> idcs;
			for (int i = 0; i < c.m_E.rows(); i++) {
				if (i == 0 && (c.m_V(c.m_E(i,0), 2) < z || !use_z)) {
					m_V_end.insert(m_V_correspondance[c_idx].at(c.m_E(i, 0)));
					m_V_end.insert(m_V_correspondance[c_idx].at(c.m_E(i, 1)));
				}
				else if (i == c.m_E.rows() - 1 && (c.m_V(c.m_E(i, 1), 2) < z || !use_z)) {
					m_V_end.insert(m_V_correspondance[c_idx].at(c.m_E(i, 0)));
					m_V_end.insert(m_V_correspondance[c_idx].at(c.m_E(i, 1)));
				}
			}
		}
	}

}

void curve_network::display(const std::string& name, bool with_pt_cloud, bool enabled) {

	//std::cout << "Displaying curve network" << std::endl;

	polyscope::registerCurveNetwork(name, m_V, m_E);
	if (with_pt_cloud) {
		polyscope::registerPointCloud(name, m_V);
		polyscope::getPointCloud(name)->setEnabled(false);
	}
	polyscope::getCurveNetwork(name)->setEnabled(enabled);
}

void curve_network::display_per_curve_with_radius(const std::string& name, const Eigen::VectorXd& rs) {
	std::cout << "displaying curve network with radiis" << std::endl;
	for (int i = 0; i < m_curves.size(); i++) {
		std::string numeral;
		numeral = i >= 10 ? "" : "0";
		auto a = polyscope::registerCurveNetwork(name + " " + numeral + std::to_string(i), m_curve_set.m_all_curves[m_curves[i]].m_V, m_curve_set.m_all_curves[m_curves[i]].m_E);
		a->setRadius(rs(i), false);
	}
}

void curve_network::display_per_curve_with_radius(const std::string& name, const std::vector<double>& rs) {
	std::cout << "displaying curve network with radiis" << std::endl;
	std::vector<double> rstest = {1,2,3,4,5};

	Eigen::VectorXd b = Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(rs.data(), rs.size());
	std::array<double, 3> random_color = { {polyscope::randomUnit(), polyscope::randomUnit(), polyscope::randomUnit()} };
	
	for (int i = 0; i < m_curves.size(); i++) {

		std::vector<std::array<double, 3>> randColor(m_curve_set.m_all_curves[m_curves[i]].m_E.size() / 2);
		for (size_t j = 0; j < m_curve_set.m_all_curves[m_curves[i]].m_E.size() / 2; j++) {
			randColor[j] = random_color;
		}

		std::string numeral;
		numeral = (i + 1) >= 10 ? "" : "0";
		auto a = polyscope::registerCurveNetwork(name + " " + numeral + std::to_string(i + 1), m_curve_set.m_all_curves[m_curves[i]].m_V, m_curve_set.m_all_curves[m_curves[i]].m_E);
		a->setRadius(b(i), false);
		a->addEdgeColorQuantity("random color", randColor);
	}
	std::cout << "displaying curve network with radiis : done" << std::endl;
}

void curve_network::scale(double s)
{
	m_V *= s;
	m_intersection_set.scale(s);
	update_curves();
}

void curve_network::generate_elastic_energy_hessian(bool do_update) {

	//std::cout << "Generating elastic energy hessian" << std::endl;

	double epsilon = 1e-8;

	bool fix_ends = true;

	if(do_update)
		update();

	// the hessian : 
	int n = m_V.rows() * 3;
	m_H = SpMat(n, n);

	// For every curve in the curve network, we will compute an hessian. 
	// The problem is transfering it into the hessian of the network.
	// To make a correspondance between vertices from a curve
	// and vertices in the network, we use m_V_correspondance

	std::vector<Eigen::Triplet<double>> triplets;
	for (int c_idx = 0; c_idx < m_curves.size(); c_idx++) {
		curve c = m_curve_set.m_all_curves[m_curves[c_idx]];

		// the hessian of curve c will be c.m_H
		c.compute_full_elastic_energy_hessian();
		
		// iterates over nonzero coefficient of sparse matrix
		for (int i = 0; i < c.m_H.outerSize(); i++) {
			for (SpMat::InnerIterator it(c.m_H, i); it; ++it) {
				int r = it.row(), c = it.col();
				double v = it.value();

				int r_flr = floor(r / 3), c_flr = floor(c / 3);
				int r_rmdr = r % 3, c_rmdr = c % 3;

				int r_corresp_idx = m_V_correspondance[c_idx].at(r_flr), c_corresp_idx = m_V_correspondance[c_idx].at(c_flr);

				int row_idx = 3 * r_corresp_idx + r_rmdr, col_idx = 3 * c_corresp_idx + c_rmdr;
			
				triplets.push_back({ row_idx, col_idx, v });
			}
		}
	}
	m_H.setFromTriplets(triplets.begin(), triplets.end());

	//std::cout << "Generating elastic energy hessian : done" << std::endl;

}

void curve_network::compute_per_curve_hessian_yes_no()
{
	//std::cout << ". compute_per_curve_hessian_yes_no()" << std::endl;

	generate_elastic_energy_hessian(true);
	int n = m_V.rows() * 3;

	for (int c_idx = 0; c_idx < m_curves.size(); c_idx++) {
		// m_IP -> current_curve_IP
		std::set<int> current_curve_IP;
		std::map<int, int> correspo = m_V_correspondance[c_idx];
		for (auto ip : m_IP) {

			for (auto e : correspo) {
				if (e.second == ip) {
					current_curve_IP.insert(e.first);
					break;
				}
			}
		}
		SpMat tmp;
		m_H_yes_stretch.push_back(tmp), m_H_yes_bend.push_back(tmp); m_H_no_stretch.push_back(tmp); m_H_no_bend.push_back(tmp); m_H_penalty.push_back(tmp);

		m_curve_set.m_all_curves[m_curves[c_idx]].compute_yes_no_stretching_hessian(m_H_yes_stretch[c_idx], m_H_no_stretch[c_idx], current_curve_IP);
		m_curve_set.m_all_curves[m_curves[c_idx]].compute_yes_no_bending_hessian(m_H_yes_bend[c_idx], m_H_no_bend[c_idx], current_curve_IP);
		m_curve_set.m_all_curves[m_curves[c_idx]].compute_penalty_hessian(m_H_penalty[c_idx]);

		// im gonna use those 2 vectors so that i can have a loop for writing all those hessians
		std::vector<SpMat> loop_small_h = { m_H_yes_stretch[c_idx], m_H_no_stretch[c_idx], m_H_yes_bend[c_idx], m_H_no_bend[c_idx], m_H_penalty[c_idx] };
		std::vector<std::string> affixes = { "yes_stretch", "no_stretch", "yes_bend", "no_bend", "penalty" };
		SpMat big_h(n, n);

		// export those hessians ... 
		for (int vec_loop = 0; vec_loop < loop_small_h.size(); vec_loop++) {
			std::vector<Eigen::Triplet<double>> triplets;
			for (int i = 0; i < loop_small_h[vec_loop].outerSize(); i++) {
				for (SpMat::InnerIterator it(loop_small_h[vec_loop], i); it; ++it) {
					int r = it.row(), c = it.col();
					double v = it.value();
					int r_flr = floor(r / 3), c_flr = floor(c / 3);
					int r_rmdr = r % 3, c_rmdr = c % 3;
					int r_corresp_idx = m_V_correspondance[c_idx].at(r_flr), c_corresp_idx = m_V_correspondance[c_idx].at(c_flr);
					int row_idx = 3 * r_corresp_idx + r_rmdr, col_idx = 3 * c_corresp_idx + c_rmdr;
					triplets.push_back({ row_idx, col_idx, v });
				}
			}
			loop_small_h[vec_loop] = SpMat(n, n);
			loop_small_h[vec_loop].setFromTriplets(triplets.begin(), triplets.end());
		}

		m_H_yes_stretch[c_idx] = loop_small_h[0];
		m_H_no_stretch[c_idx] = loop_small_h[1];
		m_H_yes_bend[c_idx] = loop_small_h[2];
		m_H_no_bend[c_idx] = loop_small_h[3];
		m_H_penalty[c_idx] = loop_small_h[4];
	}
	//std::cout << "call : compute_per_curve_hessian_yes_no : done" << std::endl;
}

void curve_network::compute_per_curve_length()
{
	m_curve_len = Eigen::VectorXd::Zero(m_curves.size());

	for (int i = 0; i < m_curves.size(); i++) {
		m_curve_len(i) = m_curve_set.m_all_curves[m_curves[i]].get_length();
	}
}

void curve_network::resample_cn(double l) {
	//std::cout << "Resampling curve network" << std::endl;

	if (l == 0)
		l = gc_resample_length;

	Eigen::MatrixXd inter_vertices;
	for (int c_idx : m_curves) {
		m_curve_set.m_all_curves[c_idx].update_inter_pts_segments(m_intersection_set);

		// gather every intersection vertex position for curve with c_idx
		int nb_inter = m_intersection_set.get_intersection_vertices(c_idx, inter_vertices);
		m_curve_set.m_all_curves[c_idx].resample(l, inter_vertices);

		m_curve_set.m_all_curves[c_idx].update_inter_pts_segments(m_intersection_set);
	}

	//std::cout << "Resampling curve network : done" << std::endl;
}

void curve_network::set_m_V_initial() {

	// need to go into each curve currently listed in m_curves and construct m_V_initial from THEIR m_V_initial + our m_V_correspondance
	m_V_initial = m_V;

	for (int i = 0; i < m_curves.size(); i++){ 
		int c_idx = m_curves[i];
		const curve& c = m_curve_set.m_all_curves[c_idx];
		const auto& mvc = m_V_correspondance[i];
		
		for (auto const& e : mvc) {
			m_V_initial.row(e.second) = c.m_V.row(e.first);
		}
	}
}

void curve_network::print_m_V_correspondance()
{
	// std::vector<std::map<int, int>> m_V_correspondance;

	std::cout << std::endl;
	std::cout << "###############################" << std::endl;
	std::cout << "Printing m_V_correspondance" << std::endl;

	std::cout << "Nb of active curves : " << m_curves.size() << std::endl;
	std::cout << "Nb of active curves according to m_V_correspondance : " << m_V_correspondance.size() << std::endl;
	std::cout<<std::endl;
	for (int c_idx = 0; c_idx < m_curves.size(); c_idx++) {
		std::cout << "curve " << m_curves[c_idx] << " (" << c_idx + 1 << "/" << m_curves.size() << ")" << std::endl;

		for (std::map<int, int>::iterator iter = m_V_correspondance[c_idx].begin(); iter != m_V_correspondance[c_idx].end(); ++iter)
		{
			int k = iter->first;
			//ignore value
			//Value v = iter->second;

			std::cout << iter->first << " -> " << iter->second << std::endl;
		}
		std::cout << std::endl;

	}

	std::cout << "Printing m_V_correspondance : done" << std::endl;
	std::cout << "###############################" << std::endl;
	std::cout << std::endl;
}

void curve_network::add_curve_with_cuts(curve c)
{
	// only called for figure 11d
	std::cout << std::endl << "ADDING CURVE " << m_curves.size() << std::endl;

	auto compare_idx = [](std::tuple<int, std::vector<Eigen::Vector3d>> a, std::tuple<int, std::vector<Eigen::Vector3d>> b) {
		return (std::get<0>(a) < std::get<0>(b));
	};

	auto remove_duplicate = [](std::vector<Eigen::Vector3d>& vec) {

		std::vector<Eigen::Vector3d> vec2;
		for (int i = 0; i < vec.size(); i++) {
			bool already_in = false;
			for (int j = 0; j < vec2.size(); j++) {
				if ((vec[i] - vec2[j]).norm() < 1e-9) {
					already_in = true;
					break;
				}
			}
			if (!already_in) {
				vec2.push_back(vec[i]);
			}
		}

		return vec2;
	};

	c.check_closed();
	
	int new_curve_idx = m_curve_set.add_curve(c) - 1;
	m_curves.push_back(new_curve_idx);

	std::vector<std::tuple<int, std::vector<Eigen::Vector3d>>> pts_to_add_per_edge; // tuples (idx edge, vector of new points on the edge)
	std::vector<Eigen::Vector3d> cut_pts;
	bool cut_done = false;
	double epsilon = 1e-8;
	int tmp_count = 0;
	for (int idx = 0; idx < m_curves.size() - 1; idx++) {
		int c_idx = m_curves[idx];

		std::vector<std::tuple<int, std::vector<Eigen::Vector3d>>> pts_to_add_per_edge_old_curve;

		// for each segment in the new curve
		for (int e1 = 0; e1 < m_curve_set.m_all_curves[m_curves.size() - 1].m_E.rows(); e1++) {

			// for each segment in the current curve of the curve network we want to compare to
			for (int e2 = 0; e2 < m_curve_set.m_all_curves[c_idx].m_E.rows(); e2++) {

				Eigen::Vector3d v0 = c.m_V.row(c.m_E(e1, 0)), v1 = c.m_V.row(c.m_E(e1, 1)), v2 = m_curve_set.m_all_curves[c_idx].m_V.row(m_curve_set.m_all_curves[c_idx].m_E(e2, 0)), v3 = m_curve_set.m_all_curves[c_idx].m_V.row(m_curve_set.m_all_curves[c_idx].m_E(e2, 1));
				bool should_intersect = false;

				if ((v1 - v0).norm() <= 1e-8 || (v3 - v2).norm() <= 1e-8) {
					continue;
				}

				if (((v1 - v0).cross(v3 - v2)).norm() <= 1e-9) {
					if (c_idx == 16 && m_curves.size() - 1 == 27) {
						tmp_count++;
						continue;
					}
				}

				if (fabs((v0 - v2).norm()) > epsilon && fabs((v0 - v3).norm()) > epsilon && fabs((v1 - v2).norm()) > epsilon && fabs((v1 - v3).norm()) > epsilon) {
					bool same_plane = check_coplanarity_lines(v0, v1, v2, v3, 0.002);
					if (c_idx == 13 && m_curves.size() == 31) {
						Eigen::Vector3d p0 = (v1 - v0).normalized();
						Eigen::Vector3d p2 = (v3 - v2).normalized();
						double val = abs((p0.cross(p2)).dot((v2 - v0).normalized()));
						if (val < 1e-3) {
							std::cout << "same PLANE ?!" << val << std::endl;
							std::cout<<"  "<< v0.transpose() << " -- " << v1.transpose() << std::endl;
							std::cout << "  " << v2.transpose() << " -- " << v3.transpose() << std::endl;

						}

						if (e2 >= m_curve_set.m_all_curves[c_idx].m_E.rows() - 2) {
							std::cout << "last 2 edges val : " << val << std::endl;
						}

					}
					if (same_plane) {
						Eigen::Vector3d inter_pt;
						if (c_idx == 13 && m_curves.size() == 31) {
							std::cout << "SAME PLANE ?!" << std::endl;
						}
						double t, s;
						find_line_intersection(v0, v1, v2, v3, inter_pt, t, s);
						if (c_idx == 13 && m_curves.size() == 31) {
							std::cout << "t " <<t<<"s "<<s<< std::endl;
						}

						bool intersection_exists = (s >= 0 && s <= 1 && t >= 0 && t <= 1);
						if (intersection_exists) {
							int inter_idx = m_intersection_set.add_intersection(inter_pt, c_idx, m_curves.size() - 1);
							inter_pt = m_intersection_set.m_all_intersections[inter_idx].m_pos;
							bool first_time_splitting = true;

							if (!(s > 1 - 1e-8 || s < 1e-8)) {
								for (int i = 0; i < pts_to_add_per_edge.size(); i++) {
									if (std::get<0>(pts_to_add_per_edge[i]) == e1) {
										std::get<1>(pts_to_add_per_edge[i]).push_back(inter_pt);
										first_time_splitting = false;
										break;
									}
								}
								if (first_time_splitting) {
									std::vector<Eigen::Vector3d> tmp;
									tmp.push_back(inter_pt);
									pts_to_add_per_edge.push_back(std::make_tuple(e1, tmp));
								}
							}

							// and on the old curve IF it's not the new curve (last curve of the m_curves vector
							first_time_splitting = true;
							if (!(t > 1 - 1e-8 || t < 1e-8)) {
								for (int i = 0; i < pts_to_add_per_edge_old_curve.size(); i++) {
									if (std::get<0>(pts_to_add_per_edge_old_curve[i]) == e2) {
										std::get<1>(pts_to_add_per_edge_old_curve[i]).push_back(inter_pt);
										first_time_splitting = false;
										break;
									}
								}
								if (first_time_splitting) {
									std::vector<Eigen::Vector3d> tmp;
									tmp.push_back(inter_pt);
									pts_to_add_per_edge_old_curve.push_back(std::make_tuple(e2, tmp));
								}
							}

							std::vector<int> curves_vec;
							curves_vec.push_back(m_curves.size() - 1);
							curves_vec.push_back(c_idx);

							if (c_idx >= 0 && c_idx <= 5) {
								cut_pts.push_back(inter_pt);
							}
						}
					}
				}

				// they intersect at an already existing vertex.
				else if (fabs((v0 - v2).norm()) <= epsilon) {
					m_intersection_set.add_intersection(v0, c_idx, m_curves.size() - 1);
					m_curve_set.m_all_curves[m_curves.size() - 1].m_inter_pts.insert(c.m_E(e1, 0));
					m_curve_set.m_all_curves[c_idx].m_inter_pts.insert(m_curve_set.m_all_curves[c_idx].m_E(e2, 0));
				}
				else if (fabs((v1 - v2).norm()) <= epsilon) {
					m_intersection_set.add_intersection(v1, c_idx, m_curves.size() - 1);
					m_curve_set.m_all_curves[m_curves.size() - 1].m_inter_pts.insert(c.m_E(e1, 1));
					m_curve_set.m_all_curves[c_idx].m_inter_pts.insert(m_curve_set.m_all_curves[c_idx].m_E(e2, 0));
				}
				else if (fabs((v0 - v3).norm()) <= epsilon) {
					m_intersection_set.add_intersection(v0, c_idx, m_curves.size() - 1);
					m_curve_set.m_all_curves[m_curves.size() - 1].m_inter_pts.insert(c.m_E(e1, 0));
					m_curve_set.m_all_curves[c_idx].m_inter_pts.insert(m_curve_set.m_all_curves[c_idx].m_E(e2, 1));
				}
				else if (fabs((v1 - v3).norm()) <= epsilon) {
					if (c_idx == 16 && m_curves.size() - 1 == 30) { std::cout << "say smth" << std::endl; }
					//std::cout << "smth is up at " << v1.transpose() << " with " << c_idx << " " << m_curves.size() - 1 << std::endl;
					m_intersection_set.add_intersection(v1, c_idx, m_curves.size() - 1);
					//std::cout << "2 vertices kinda close ?" << std::endl;
					m_curve_set.m_all_curves[m_curves.size() - 1].m_inter_pts.insert(c.m_E(e1, 1));
					m_curve_set.m_all_curves[c_idx].m_inter_pts.insert(m_curve_set.m_all_curves[c_idx].m_E(e2, 1));
				}
			}
		}
		// would use a lambda, but m_curves is not a variable ...
		
		std::sort(pts_to_add_per_edge_old_curve.begin(), pts_to_add_per_edge_old_curve.end(), compare_idx);

		for (int i = pts_to_add_per_edge_old_curve.size() - 1; i >= 0; i--) {
			std::tuple<int, std::vector<Eigen::Vector3d>> tupl = pts_to_add_per_edge_old_curve[i];
			int e = std::get<0>(tupl);
			std::vector<Eigen::Vector3d> l = std::get<1>(tupl);
			l = remove_duplicate(l);
			m_curve_set.m_all_curves[c_idx].split_edge(e, l);
		}	std::cout << std::endl;
	}

	std::sort(pts_to_add_per_edge.begin(), pts_to_add_per_edge.end(), compare_idx);

	for (int i = pts_to_add_per_edge.size() - 1; i >= 0; i--) {
		std::tuple<int, std::vector<Eigen::Vector3d>> tupl = pts_to_add_per_edge[i];
		int e = std::get<0>(tupl);
		std::vector<Eigen::Vector3d> l = std::get<1>(tupl);
		l = remove_duplicate(l);
		m_curve_set.m_all_curves[m_curves.size() - 1].split_edge(e, l);
	}

	for (Eigen::Vector3d p : cut_pts) {
		m_curve_set.m_all_curves[m_curves.size() - 1].cut(p);
	}
}

void curve_network::write_m_V(const std::string& path) {
	std::ofstream f;
	f.open(path);

	for (int v_idx = 0; v_idx < m_V.rows(); v_idx++) {
		f << to_string_with_precision(m_V.row(v_idx), 15) << "\n";
	}
}

void curve_network::read_m_V(const std::string& path) {
	std::cout << std::endl << "Reading m_V from " << path << std::endl;
	
	Eigen::MatrixXd new_m_V(m_V.rows(), m_V.cols());

	std::ifstream file(path);
	if (!file) {
		std::cerr << "Failed to load vertex data\n";
		return;
	}

	std::string s;
	std::string subs;
	int i = 0;

	while (std::getline(file, s)) {
		std::stringstream ss;
		ss.str(s);
		bool new_curve_line = false;


		while (ss >> subs) {
			// taking 1 of 3 coordinates at a time
			new_m_V(floor(i/3), i % 3) = atof(subs.c_str());
			i += 1;
		}
	}
	
	m_V = new_m_V;
	update_curves();

	return;
}

void curve_network::write_curves(const std::string& path, std::vector<double> a)
{
	std::ofstream f, f2;
	f.open("output/greedy/" + path + "_full_curves.txt");
	f2.open("output/greedy/" + path + "_optim_curves.txt");
	//    myfile << "Writing this to a file.\n";
	int counter = 0;


	for (curve c : m_curve_set.m_all_curves) {
		f << "CURVE " << counter << "\n";
		for (int i = 0; i < c.m_E.rows(); i++) {
			f << to_string_with_precision(c.m_V.row(c.m_E(i, 0)), 15) << "\n";
		}
		f << to_string_with_precision(c.m_V.row(c.m_E(c.m_E.rows() - 1, 1)), 15) << "\n";
	}

	for (curve c : m_curve_set.m_all_curves) {
		
		if (a[counter] > 0) {
			f2 << "CURVE " << counter << "\n";
			for (int i = 0; i < c.m_E.rows(); i++) {
				f2 << to_string_with_precision(c.m_V.row(c.m_E(i, 0)), 15) << "\n";
			}
			f2 << to_string_with_precision(c.m_V.row(c.m_E(c.m_E.rows() - 1, 1)), 15) << "\n";
		}
		counter++;
	}
	f.close();
	f2.close();

}

void curve_network::write_curves(const std::string& path, Eigen::VectorXd a)
{
	std::ofstream f, f2;
	f.open("output/greedy/" + path + "_full_curves.txt");
	f2.open("output/greedy/" + path + "_optim_curves.txt");
	//    myfile << "Writing this to a file.\n";
	int counter = 0;


	for (curve c : m_curve_set.m_all_curves) {
		f << "CURVE " << counter << "\n";
		for (int i = 0; i < c.m_E.rows(); i++) {
			f << to_string_with_precision(c.m_V.row(c.m_E(i, 0)), 15) << "\n";
		}
		f << to_string_with_precision(c.m_V.row(c.m_E(c.m_E.rows() - 1, 1)), 15) << "\n";
	}

	for (curve c : m_curve_set.m_all_curves) {

		if (a(counter) > 0) {
			f2 << "CURVE " << counter << "\n";
			for (int i = 0; i < c.m_E.rows(); i++) {
				f2 << to_string_with_precision(c.m_V.row(c.m_E(i, 0)), 15) << "\n";
			}
			f2 << to_string_with_precision(c.m_V.row(c.m_E(c.m_E.rows() - 1, 1)), 15) << "\n";
		}
		counter++;
	}
	f.close();
	f2.close();

}

void curve_network::write_curves(const std::string& path)
{
	std::ofstream f, f2;
	f.open(path);
	int counter = 0;

	for (int i = 0; i < m_curves.size(); i++) {
		curve c = m_curve_set.m_all_curves[m_curves[i]];
		f << "CURVE " << i << "\n";
		for (int i = 0; i < c.m_E.rows(); i++) {
			f << to_string_with_precision(c.m_V.row(c.m_E(i, 0)), 15) << "\n";
		}
		f << to_string_with_precision(c.m_V.row(c.m_E(c.m_E.rows() - 1, 1)), 15) << "\n";
	}

	f.close();
}

void curve_network::lock_initial_cn() {
	for (int i = 0; i < m_curves.size(); i++) {
		m_curve_set.m_all_curves[m_curves[i]].compute_init_lengths();
		m_curve_set.m_all_curves[m_curves[i]].set_m_V_initial();
	}
}
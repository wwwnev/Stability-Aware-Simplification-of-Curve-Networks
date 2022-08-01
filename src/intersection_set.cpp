#include <iostream>

#include "intersection_set.h"
#include "intersection.h"

int intersection_set::add_intersection(const Eigen::Vector3d& pos, int c_id1, int c_id2)
{
	// check if there's already an intersection at pos

	for (int i = 0; i < m_all_intersections.size(); i++) {
		if ((m_all_intersections[i].m_pos - pos).norm() < 5e-3) {
			
			m_all_intersections[i].add_curve(c_id1);
			if (c_id2 != -1)
				m_all_intersections[i].add_curve(c_id2);
			return i;
		}
	} 

	// else create an intersection
	if (c_id2 != -1)
		m_all_intersections.push_back(intersection(pos, c_id1, c_id2));
	return m_all_intersections.size() - 1;
}

int intersection_set::remove_curve(int c_id)
{
	for (auto it = m_all_intersections.begin(); it != m_all_intersections.end(); ) {
		int nb_inter_curves = it->remove_curve(c_id);
		if (nb_inter_curves <= 1)
			it = m_all_intersections.erase(it);
		else
			++it;
	}
	return m_all_intersections.size();
}

int intersection_set::remove_intersection(const Eigen::Vector3d& pos)
{
	for (auto it = m_all_intersections.begin(); it != m_all_intersections.end(); ) {
		if ((it->m_pos - pos).norm() < 1e8) {
			it = m_all_intersections.erase(it);
		}
		else
			it++;
	}

	return 0;
}

int intersection_set::get_intersection_vertices(int c_idx, Eigen::MatrixXd& intersection_vertices)
{
	intersection_vertices = Eigen::MatrixXd(0, 3);
	int inter_count = 0;
	for (intersection inter : m_all_intersections) {
		const bool is_in = inter.m_curve_ids.find(c_idx) != inter.m_curve_ids.end();
		if (inter.m_curve_ids.find(c_idx) != inter.m_curve_ids.end()) {
			inter_count++;
			intersection_vertices.conservativeResize(intersection_vertices.rows() + 1, 3);
			intersection_vertices.row(intersection_vertices.rows() - 1) = inter.m_pos;
		}
	}

	return inter_count;
}

int intersection_set::get_nb_intersections(int c_idx)
{
	int nb_inter = 0;
	for (intersection inter : m_all_intersections) {
		if (inter.m_curve_ids.find(c_idx) != inter.m_curve_ids.end())
			nb_inter++;
	}
	return nb_inter;
}

Eigen::VectorXd intersection_set::merge_3_intersections(std::vector<int> intr_ids)
{
	Eigen::Vector3d new_pos = Eigen::Vector3d::Zero();
	std::set<int> crvs;
	std::sort(intr_ids.begin(), intr_ids.end());

	for (int i = intr_ids.size() - 1; i >= 0; i--) {
		int i_idx = intr_ids[i];
		new_pos += m_all_intersections[i_idx].m_pos;
		for (auto a : m_all_intersections[i_idx].m_curve_ids)
			crvs.insert(a);
		m_all_intersections.erase(m_all_intersections.begin() + i_idx);
	}
	new_pos = new_pos / 3;
	std::vector<int> crvs2;
	for (auto a : crvs)
		crvs2.push_back(a);

	intersection new_intersection(new_pos, crvs2[0], crvs2[1]);
	new_intersection.add_curve(crvs2[2]);
	m_all_intersections.push_back(new_intersection);
	return new_pos;
}

void intersection_set::print_intersection_set() const
{
	std::cout << "intersection set :" << std::endl;
	int count = 0;
	for (intersection inter : m_all_intersections) {
		std::cout << "inter " << count<<" : ";
		for (int c_id : inter.m_curve_ids)
			std::cout << c_id << ", ";
		std::cout << "("<<inter.m_pos.transpose()<<")"<<std::endl;
		count++;
	}
}

void intersection_set::scale(double s)
{
	for (int i = 0; i < m_all_intersections.size(); i++) {
		m_all_intersections[i].m_pos *= s;
	}
}

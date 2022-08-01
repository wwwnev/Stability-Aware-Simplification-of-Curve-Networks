#pragma once
#define _USE_MATH_DEFINES
#include <vector>

#include <Eigen/Core>

#include "intersection.h"

class intersection_set
{
public:
	std::vector<intersection> m_all_intersections;

	int add_intersection(const Eigen::Vector3d& pos, int c_id1, int c_id2 = -1);

	int remove_curve(int c_id);

	int remove_intersection(const Eigen::Vector3d& pos);

	int get_intersection_vertices(int c_idx, Eigen::MatrixXd& intersection_vertices);

	int get_nb_intersections(int c_idx);

	Eigen::VectorXd merge_3_intersections(std::vector<int> intr_ids);

	void print_intersection_set() const;
	void scale(double s);

};
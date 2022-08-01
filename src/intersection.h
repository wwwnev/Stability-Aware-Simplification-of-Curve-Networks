#pragma once
//#define _USE_MATH_DEFINES
//#include <math.h>
//#include <string>
//#include <cstdlib>
//#include <cmath>
#include <vector>
#include <set>

#include <Eigen/Core>

class intersection
{
public:
	Eigen::Vector3d m_pos;
	std::set<int> m_curve_ids;

	intersection(const Eigen::Vector3d& pos, int c_id1, int c_id2);

	int add_curve(int c_id);
	int remove_curve(int c_id); // returns the new size of m_curve_ids
};


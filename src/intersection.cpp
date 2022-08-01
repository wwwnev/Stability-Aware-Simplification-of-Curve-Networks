#include "curve.h"
#include "aux_tools.h"
#include "intersection.h"

intersection::intersection(const Eigen::Vector3d& pos, int c_id1, int c_id2)
{
	m_pos = pos;
	m_curve_ids.insert(c_id1);
	m_curve_ids.insert(c_id2);
}

int intersection::add_curve(int c_id)
{
	m_curve_ids.insert(c_id);
	return m_curve_ids.size();
}

int intersection::remove_curve(int c_id)
{
	for (auto it = m_curve_ids.begin(); it != m_curve_ids.end(); ) {
		if (*it == c_id)
			it = m_curve_ids.erase(it);
		else
			++it;
	}


	return m_curve_ids.size();
}

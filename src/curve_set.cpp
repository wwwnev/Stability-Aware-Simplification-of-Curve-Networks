#include "curve_set.h"

int curve_set::add_curve(curve c)
{
	m_all_curves.push_back(c);
	return m_all_curves.size();
}

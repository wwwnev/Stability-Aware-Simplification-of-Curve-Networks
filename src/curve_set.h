#pragma once
#define _USE_MATH_DEFINES
#include <vector>

#include <Eigen/Core>

#include "curve.h"

class curve_set
{
public:
	std::vector<curve> m_all_curves;
	int add_curve(curve c);
};

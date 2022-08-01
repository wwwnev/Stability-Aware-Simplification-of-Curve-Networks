#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include "cmath"
#include "iostream"
#include "vector"

#include <Eigen/Core>

#include "curve_network.h"
#include "surface.h"

void generate_stadium(curve_network& cn, surface& s, int n = 1);
void generate_hill(curve_network& cn, surface& s, int n = 3);
void generate_bunny(curve_network& cn, surface& s, int n);
void generate_tent(curve_network& cn, surface& s, int n);
void generate_your_example(curve_network& cn, surface& s, int n);
#pragma once

#include <string>

#include "curve_network.h"
#include "surface.h"
void load_curves_for_display(const std::string& path, curve_network& cn);

void load_curves_elim_duplicate(const std::string& path, curve_network& cn, double epsilon = 1e-5, double scale = 1.0);

void load_curves(const std::string& path, curve_network& cn, double scale = 1.0);
void load_curves(const std::string& path, curve_network& cn, surface& surf, double scale = 1.0);

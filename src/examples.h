#pragma once
#include "curve_network.h"

void roof(curve_network& cn);

void stadium(curve_network& cn, int nb_curves = 10);

void shell(curve_network& cn);

void arcshell(curve_network& cn);

void kagome(curve_network& cn);

void bunny(curve_network& cn, int nb_curves = 3);

void hill(curve_network& cn, int nb_curves = 3);

void tent(curve_network& cn, int nb_curves);

void tower(curve_network& cn);

void your_example(curve_network& cn, std::string name, std::string in_path, std::string out_path, int nb_curves, double budget);

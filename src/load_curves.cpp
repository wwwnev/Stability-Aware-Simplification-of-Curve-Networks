#include <sstream>
#include <vector>
#include <stdlib.h>
#include <Eigen/Core>
#include <iostream>
#include <fstream>

#include "curve.h"
#include "load_curves.h"

// https://codereview.stackexchange.com/questions/189536/reading-vertex-data-from-a-file
void load_curves_for_display(const std::string& path, curve_network& cn) {
    std::cout << std::endl << "Loading curves from " << path << std::endl;

    std::vector<curve> curves;

    std::ifstream file(path);
    if (!file) {
        std::cerr << "Failed to load vertex data\n";
        return;
    }

    int curve_idx = -1;
    std::string s;
    std::string subs;
    int i = 0;
    Eigen::MatrixXd v0(3, 1), v1(3, 1);
    while (std::getline(file, s)) {
        std::stringstream ss;
        ss.str(s);
        bool new_curve_line = false;


        while (ss >> subs) {
            if (subs == "CURVE") {
                curves.push_back(curve());
                curve_idx += 1;
                i = 0;
                break;
            }
            else if (curve_idx >= 0) {
                // taking 1 of 3 coordinates at a time
                if (i % 3 == 0) {
                    v1(0) = atof(subs.c_str());
                }
                else if (i % 3 == 1) {
                    v1(1) = atof(subs.c_str());
                }
                // when it's the third coordinate, add an edge
                else if (i % 3 == 2) {
                    v1(2) = atof(subs.c_str());
                    // except if it's the first vertex, need a second one to form an edge !
                    if (i > 2) {
                        curves[curve_idx].add_edge(v0, v1);
                    }
                    v0 = v1;
                }
                i += 1;
            }
        }
    }

    // filling the curve network cn
    for (curve c : curves) {
        //c.m_V *= scale;
        cn.add_curve_for_display(c);
    }


    return;
}

void load_curves_elim_duplicate(const std::string& path, curve_network& cn, double epsilon, double scale)
{
    std::cout << std::endl << "Loading curves from " << path << std::endl;

    std::vector<curve> curves;

    std::ifstream file(path);
    if (!file) {
        std::cerr << "Failed to load vertex data\n";
        return;
    }

    int curve_idx = -1;
    std::string s;
    std::string subs;
    int i = 0;
    //double x0 = 0.0, y0 = 0.0, z0 = 0.0, x1 = 0.0, y1 = 0.0, z1 = 0.0;
    Eigen::MatrixXd v0(3, 1), v1(3, 1);
    while (std::getline(file, s)) {
        //std::cout << s << std::endl;
        std::stringstream ss;
        ss.str(s);
        bool new_curve_line = false;


        while (ss >> subs) {
            if (subs == "CURVE") {
                curves.push_back(curve());
                curve_idx += 1;
                i = 0;
                break;
            }
            else if (curve_idx >= 0) {
                // taking 1 of 3 coordinates at a time
                if (i % 3 == 0) {
                    v1(0) = atof(subs.c_str());
                }
                else if (i % 3 == 1) {
                    v1(1) = atof(subs.c_str());
                }
                // when it's the third coordinate, add an edge
                else if (i % 3 == 2) {
                    v1(2) = atof(subs.c_str());
                    // except if it's the first vertex, need a second one to form an edge !
                    if (i <= 2) {
                        v0 = v1;
                    }
                    else if ((v1 - v0).norm() > epsilon) {
                        curves[curve_idx].add_edge(v0, v1);
                        v0 = v1;
                    }
                }
                i += 1;
            }
        }
    }

    // filling the curve network cn
    i = 0;
    for (curve c : curves) {
        c.m_V *= scale;
        std::cout << "c " << i << ", " << c.m_E.rows() << ", " << c.m_V.rows() << std::endl;
        cn.add_curve(c);
    }

    return;
}

void load_curves(const std::string& path, curve_network& cn, double scale) {
	std::cout <<std::endl<< "Loading curves from " << path<< std::endl;

    std::vector<curve> curves;

    std::ifstream file(path);
    if (!file) {
        std::cerr << "Failed to load vertex data\n";
        return;
    }

    int curve_idx = -1;
    std::string s;
    std::string subs;
    int i = 0;
    //double x0 = 0.0, y0 = 0.0, z0 = 0.0, x1 = 0.0, y1 = 0.0, z1 = 0.0;
    Eigen::MatrixXd v0(3, 1), v1(3, 1);
    while (std::getline(file, s)) {
        //std::cout << s << std::endl;
        std::stringstream ss;
        ss.str(s);
        bool new_curve_line = false;
        
        
        while (ss >> subs) {
            if (subs == "CURVE") {
                curves.push_back(curve());
                curve_idx += 1;
                i = 0;
                break;
            }
            else if (curve_idx >= 0) {
                // taking 1 of 3 coordinates at a time
                if (i % 3 == 0) {
                    v1(0) = atof(subs.c_str());
                }
                else if (i % 3 == 1) {
                    v1(1) = atof(subs.c_str());
                }
                // when it's the third coordinate, add an edge
                else if (i % 3 == 2) {
                    v1(2) = atof(subs.c_str());
                    // except if it's the first vertex, need a second one to form an edge !
                    if (i > 2) {
                        curves[curve_idx].add_edge(v0,v1);
                    }
                    //std::cout << curve_idx << " V0 : " << v0.transpose() << std::endl;
                    //std::cout << curve_idx<<" V1 : " << v1.transpose() << std::endl;
                    v0 = v1;
                }
                i += 1;
            }
        }
    }

    // filling the curve network cn
    for (curve c : curves) {
        //c.m_V *= 2;
        c.m_V *= scale;
        cn.add_curve(c);
    }

    //cn.update();
    //cn.resample_cn(1);
    //cn.update();
    //cn.set_m_V_initial();

	return;
}

void load_curves(const std::string& path, curve_network& cn, surface& surf, double scale) {
    std::cout << std::endl << "Loading curves from " << path << std::endl;

    std::vector<curve> curves;

    std::ifstream file(path);
    if (!file) {
        std::cerr << "Failed to load vertex data\n";
        return;
    }

    int curve_idx = -1;
    std::string s;
    std::string subs;
    int i = 0;
    //double x0 = 0.0, y0 = 0.0, z0 = 0.0, x1 = 0.0, y1 = 0.0, z1 = 0.0;
    Eigen::MatrixXd v0(3, 1), v1(3, 1);
    while (std::getline(file, s)) {
        //std::cout << s << std::endl;
        std::stringstream ss;
        ss.str(s);
        bool new_curve_line = false;


        while (ss >> subs) {
            if (subs == "CURVE") {
                curves.push_back(curve());
                curve_idx += 1;
                i = 0;
                break;
            }
            else if (curve_idx >= 0) {
                // taking 1 of 3 coordinates at a time
                if (i % 3 == 0) {
                    v1(0) = atof(subs.c_str());
                }
                else if (i % 3 == 1) {
                    v1(1) = atof(subs.c_str());
                }
                // when it's the third coordinate, add an edge
                else if (i % 3 == 2) {
                    v1(2) = atof(subs.c_str());
                    // except if it's the first vertex, need a second one to form an edge !
                    if (i > 2) {
                        curves[curve_idx].add_edge(v0, v1);
                    }
                    v0 = v1;
                }
                i += 1;
            }
        }
    }

    for (int i = 0; i < curves.size(); i++) {
        curves[i].remove_0_edges();
        curves[i].m_V *= scale;
    }

    if (path == "input/boring_bean_100_gridshell_for_arcs.txt") {
        for (int i = 0; i < curves.size(); i++) {
            cn.add_curve_with_cuts(curves[i]); // instead of calling cn.add_curve(curves[i])
        }
        return;
    }

    int n = curves.size();
    //int n = 3;
    for (int i = 0; i < n; i++) {
        cn.add_curve(curves[i]);
    }

    return;
}
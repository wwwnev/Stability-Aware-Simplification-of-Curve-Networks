#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <cstdlib>
#include <cmath>

#include <Eigen/Core>

#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>
#include <igl/edges.h>
#include <igl/boundary_loop.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include "polyscope/point_cloud.h"


#include "curve.h"
#include "aux_tools.h"
#include "curve_network.h"

#pragma once
class surface
{
	public:
		Eigen::MatrixXd m_V;
		Eigen::MatrixXi m_F;
		Eigen::MatrixXi m_E;
		Eigen::MatrixXi m_EA;
		curve_network m_CN;
		std::string m_path;
		std::vector<curve> m_geodesics;

		surface(std::string path, bool rewrite);
		surface(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

		void generate_edge_adjacency_matrix();
		int edge_index(int start, int end);
		void add_random_geodesic_curve(double deadzone_angle = 0.0);
		void generate_random_geodesic_curve(curve& ccc, double deadzone_angle = 0.0);

		void generate_geodesic_curve(curve& ccc, double t = 0, double angle = M_PI / 2, int boundary_choice = 0);
		void generate_geodesic_curve_from_pt(curve& ccc, int start_v_idx, double angle_t, int curve_switch = 1);
		void generate_geodesic_curve(curve& ccc, int start_idx, double angle = M_PI / 2, int boundary_choice = 0);

		void generate_horizontal_geodesic_curve_pt(curve& ccc, int v1, int v2, double& end_t, double& end_angle);
		void fix_duplicate_vertices();

		void fix_curve_almost_duplicate_vertices_on_mesh_vertices(curve& c);
		void trace_geodesic_between_2_pts(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, curve& c);
};


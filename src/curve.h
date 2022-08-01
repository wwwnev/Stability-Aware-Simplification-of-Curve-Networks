#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "set"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <polyscope/curve_network.h>

#include "constants.h"
#include "aux_tools.h"
#include "intersection_set.h"

#include "generated/computeStraightBendingEnergy.h"
#include "generated/computeEdgeStretchingEnergy.h"

#pragma once

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class curve {
	public:
		Eigen::MatrixXd m_V;
		Eigen::MatrixXd m_V_initial;
		Eigen::MatrixXi m_E;
		std::set<int> m_inter_pts;
		std::vector<std::vector<int>> m_segments;

		SpMat m_H;
		Eigen::MatrixXd m_G;

		double m_stiffness;
		double m_stiffness_modifier;

		Eigen::MatrixXd m_rot_mat;

		bool m_closed = false;
		
		std::vector<double> m_init_e_len;

		curve();

		void compute_segments();
		std::pair<int,int> add_edge(Eigen::MatrixXd v0, Eigen::MatrixXd v1);
		std::pair<int,int> add_edge(Eigen::MatrixXd v0, Eigen::MatrixXd v1, int e_precedent);
		bool split_edge(int e, std::vector<Eigen::Vector3d> inter_pts);
		bool split_edge(int e, std::set<Eigen::Vector3d> inter_pts);
		bool split_edge(int e, Eigen::Vector3d inter_pt);
		bool cut(Eigen::Vector3d p);
		double get_length();
		double compute_init_lengths();
		void build_vertex_adjacency_matrix(std::vector<std::vector<int>>& v_adj);

		void compute_full_elastic_energy_hessian();

		void compute_gradient();
		void compute_bending_gradient(Eigen::MatrixXd& grad_bend);
		void compute_stretching_gradient(Eigen::MatrixXd& grad_stretch);

		double compute_bending_energy();
		double compute_stretching_energy();

		void compute_yes_no_bending_hessian(SpMat& yes, SpMat& no, const std::set<int>& ip);
		void compute_yes_no_stretching_hessian(SpMat& yes, SpMat& no, const std::set<int>& ip);
		void compute_penalty_hessian(SpMat& penalty);

		void set_m_V_initial();
		void remove_0_edges();
		void resample(double l); // should take list of intersections as arg
		void resample(double l, const Eigen::MatrixXd& intersection_vertices);

		void replace_edge_by_vertex(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& pnew);

		void update_inter_pts_segments(const intersection_set& is);
		void close_curve();
		bool check_closed();

		bool join(curve& c);
		void replace_closest_vertex(Eigen::Vector3d new_pos);
		void smooth_new_segments(const Eigen::Vector3d& removed_vertex, const Eigen::MatrixXd& intersection_vertices);

	private:
		int add_vertex(Eigen::MatrixXd v);
};
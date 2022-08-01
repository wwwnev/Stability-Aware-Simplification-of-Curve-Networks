#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <set>

#include<Eigen/Core>
#include <Eigen/Sparse>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h> // i should try sparse sym solve next
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/Util/SelectionRule.h>

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include "constants.h"
#include "curve.h"
#include "aux_tools.h"
#include "curve_set.h"
#include "intersection_set.h"

#include "generated/computeStraightBendingEnergy.h"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class curve_network
{
	public:
		curve_set m_curve_set; // complete set of curves for the curve network.
		std::vector<int> m_curves; // subset of m_all_curves (int -> id in m_all_curves). only currently chosen curves.

		intersection_set m_intersection_set; // CURRENT intersections. depends on m_curves.

		Eigen::MatrixXd m_V; // all vertices
		Eigen::MatrixXd m_V_initial;
		Eigen::MatrixXi m_E; // all edges

		std::set<int> m_IP; // all intersection vertices. Note : contains curve net vertex indices and not singular curve vertex indices. Get those with m_V_correspondance I think...
		
		std::vector<std::map<int, int>> m_V_correspondance; // vector : index i of vector is the vertex correspondance for curve i. map iterator .first = index in curve, map iterator .second = index in CN
		std::set<int> m_V_end;
		SpMat m_H;
		std::vector<SpMat> m_H_yes_stretch;
		std::vector<SpMat> m_H_no_stretch;
		std::vector<SpMat> m_H_yes_bend;
		std::vector<SpMat> m_H_no_bend;
		std::vector<SpMat> m_H_penalty;
		double m_stiffness;
		Eigen::VectorXd m_curve_len;

		curve_network();
		void update_fast();
		void update_for_display();
		void update();
		void update_curves();
		void fix_m_V_end_bean_w_arcs(double z = 1e-5, bool use_z = false);
		void update_IP();
		void add_curve(curve c);
		void add_curve_for_display(curve c);
		std::tuple<int,int> add_edge(Eigen::MatrixXd v0, Eigen::MatrixXd v1);
		int add_vertex(Eigen::MatrixXd v);
		void remove_curve(int c_idx);
		void display(const std::string& name, bool with_pt_cloud = true, bool enabled = true);
		void display_per_curve_with_radius(const std::string&, const Eigen::VectorXd& rs);
		void display_per_curve_with_radius(const std::string&, const std::vector<double>& rs);
		
		void scale(double s);

		void generate_elastic_energy_hessian(bool do_update = true);

		double compute_energy();

		void compute_per_curve_hessian_yes_no();
		void save_per_curve_hessian_yes_no(const std::string& dir_name, const std::string& prefix = "hessian_c");
		void save_per_curve_hessian(const std::string& dir_name, const std::string& prefix = "hessian_c");

		void compute_per_curve_length();
		void save_per_curve_length(const std::string& dir_name, const std::string& file_name);

		void build_vertex_adjacency_matrix(std::vector<std::vector<int>>& v_adj);

		double hessian_smallest_eigenvector(Eigen::MatrixXd& vec_field);

		double hessian_smallest_eigenvalue();

		void resample_cn(double l = 0);

		void set_m_V_initial();

		void m_E_per_curve_wr_cn(std::vector<std::vector<int>>& m_E_per_curve);
		void save_m_E_per_curve_wr_cn(const std::string& folder_name);

		void print_m_V_correspondance();

		void merge_close_intersections();

		void print_per_curve_lambda();
		void export_hessians_and_stuff(std::string folder_name);

		void add_curve_with_cuts(curve c);

		void write_m_V(const std::string& path);
		void read_m_V(const std::string& path);

		void replace_m_V_from_file_closest(const std::string& path);

		void write_curves(const std::string& path, std::vector<double> a);
		void write_curves(const std::string& path, Eigen::VectorXd a);
		void write_curves(const std::string& path);

		void rationalize(double radius, std::vector<Eigen::Vector3d>& removed_vertices);

		void compute_gravity_force(const Eigen::VectorXd& a, Eigen::VectorXd& g);
		void lock_initial_cn();
	private:
		
		
		
};

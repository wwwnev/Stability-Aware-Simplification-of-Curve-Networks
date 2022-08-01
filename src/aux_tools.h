#pragma once

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/SparseQR>
#include "polyscope/point_cloud.h"

#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h> // i should try sparse sym solve next
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/MatOp/DenseSymShiftSolve.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/Util/SelectionRule.h>

#include <unsupported/Eigen/ArpackSupport>

typedef Eigen::SparseMatrix<double> SparseMat;
typedef Eigen::SimplicialLDLT<SparseMat> SparseChol;
typedef Eigen::ArpackGeneralizedSelfAdjointEigenSolver <SparseMat, SparseChol> Arpack;

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

//void compute_eigen_

int mod(int a, int b);
double mod(double a, double b);

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);

bool check_coplanarity_lines(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, double epsilon = 1e-6);
void find_line_intersection(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, Eigen::Vector3d& np, double& t, double&s);


void vector_to_mat(const Eigen::VectorXd& v, int rows, int cols, Eigen::MatrixXd& m);
void mat_to_vector(const Eigen::MatrixXd& m, Eigen::VectorXd& v);

bool find_line_intersection(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, Eigen::Vector3d& np);
void mat_to_matlab(const std::string& file_name, const Eigen::MatrixXd& m);
void sparse_mat_to_matlab(const std::string& file_name, const std::string& mat_name, SparseMat& m);
void std_vector_to_matlab(const std::string& path, const std::vector<int>& v);

void import_vector(const std::string& path, std::vector<double>& v);
void print_sparse_mat(const SparseMat& mat, const std::string& name = "name");

int rank_mat(const Eigen::MatrixXd& m);
int rank_mat(const SparseMat& m);
void Q_basis_mat(const Eigen::MatrixXd& m, Eigen::MatrixXd& q);
void Q_basis_mat(const SparseMat& m, SparseMat& q);

double compute_smallest_eigen_value(const SparseMat& m, int n = 1);

void convert_int_to_alphas(long long int k, int len_alphas, Eigen::VectorXd& alphas);

void reflection_pt_to_plane(Eigen::Vector3d& p, double A, double B, double C, double D);
double distance_pt_to_plane(const Eigen::Vector3d& p, double A, double B, double C, double D);
void normal_to_plane(double A, double B, double C, double D, Eigen::Vector3d& n, bool unit = true);

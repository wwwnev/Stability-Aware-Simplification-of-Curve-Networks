#include "aux_tools.h"


int mod(int a, int b) {
    return (a % b + b) % b;
}

//remove row from eigen matrix
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

    matrix.conservativeResize(numRows, numCols);
}

bool check_coplanarity_lines(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, double epsilon)
{
    Eigen::Vector3d v0 = (p1 - p0).normalized();
    Eigen::Vector3d v2 = (p3 - p2).normalized();
    double val = (v0.cross(v2)).dot((p2 - p0).normalized());
    return (val >= 0 - epsilon) && (val <= 0 + epsilon);
}

bool find_line_intersection(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, Eigen::Vector3d& np) {

    // https://mathworld.wolfram.com/Line-LineIntersection.html

    Eigen::Vector3d a = p1 - p0, b = p3 - p2, c = p2 - p0;

    // check for parallel lines
    if (a.cross(b).norm() == 0) {
        //std::cout << "Cannot compute intersection of parallel lines" << std::endl;
        return false;
    }

    double s = (c.cross(b)).dot(a.cross(b)) / (a.cross(b)).squaredNorm();

    np = p0 + a * s; 

    // if the intersection lies between p0 and p1, return true. if the intersection is outside: false.
    double epsilon = 1e-5;
    if (s>=0.0-epsilon && s<=1.0+epsilon /*&& t0>=-epsilon && t0<=1.0+epsilon*/) {
        return true;
    }
    else {
        return false;
    }
}

// .txt must be part of the name.
void mat_to_matlab(const std::string& file_name, const Eigen::MatrixXd& m)
{
    std::ofstream MyFile(file_name);
    for (int i = 0; i < m.rows(); i++) {
        for (int j = 0; j < m.cols(); j++) {
            MyFile << to_string_with_precision(m(i, j), 15);
            if (j != m.cols() - 1)
                MyFile << ",";
        }
        if (i != m.rows() - 1)
            MyFile << "\n";
    }
    MyFile.close();
}

void sparse_mat_to_matlab(const std::string& file_name, const std::string& mat_name, SparseMat& m)
{
    std::cout << "fn " << file_name << std::endl;
    std::ofstream MyFile(file_name);

    for (int k = 0; k < m.outerSize(); ++k) {
        for (SparseMat::InnerIterator it(m, k); it; ++it) {
            auto v = it.value();
            if (abs(v) < 1e-9) 
                continue;
            auto r = it.row();   // row index
            auto c = it.col();   // col index (here it is equal to k)
            MyFile << std::to_string(r + 1) + " " + std::to_string(c + 1) + " " + to_string_with_precision(v, 15) + "\n";
        }
    }
    MyFile << std::to_string(m.rows()) + " " + std::to_string(m.cols()) + " 0.0";
    MyFile.close();
}

void std_vector_to_matlab(const std::string& path, const std::vector<int>& v)
{
    // i want the path from where we are running the program, smth like : "output/folder/file_name"
    std::ofstream MyFile(path);

    for (int i = 0; i < v.size(); i++) {
        MyFile << std::to_string(v[i] + 1); // that is so dumb, it only works for exporting edge sequences to matlab, cause you need to do + 1 to indices .....
        if (i < v.size() - 1) MyFile << ",";
    }
    MyFile.close();
}

void import_vector(const std::string& path, std::vector<double>& v)
{
    // importing vector, the file should be a one liner of this format : 
    // a,b,c,d,e,f,g,h, ....... \n

    std::ifstream file(path);
    if (!file) {
        std::cerr << "Failed to load vector data\n";
        return;
    }

    std::string s;
    std::string subs;
    while (std::getline(file, s)) {
        std::stringstream ss;
        ss.str(s);

        std::string delimiter = ",";
        size_t last = 0;
        size_t next = 0;
        std::string token;
        while ((next = s.find(delimiter, last)) != std::string::npos) { 
            v.push_back(atof(s.substr(last, next - last).c_str()));
            last = next + 1; 
        }
        v.push_back(atof(s.substr(last).c_str()));
    }
    return;
}



void find_line_intersection(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, Eigen::Vector3d& np, double& t, double& s) {

    // https://mathworld.wolfram.com/Line-LineIntersection.html

    Eigen::Vector3d a = p1 - p0, b = p3 - p2, c = p2 - p0;

    // check for parallel lines
    if (a.cross(b).norm() == 0) {
        //std::cout << "Cannot compute intersection of parallel lines" << std::endl;
        return;
    }
    s = (c.cross(b)).dot(a.cross(b)) / (a.cross(b)).squaredNorm();

    np = p0 + a * s;

    // if the intersection lies between p0 and p1, return true. if the intersection is outside: false.
    double epsilon = 1e-4;
    t = (np(0) - p2(0)) / b(0); //, t1 = (np(1) - p2(1)) / b(1), t2 = (np(2) - p2(2)) / b(2);
}

void vector_to_mat(const Eigen::VectorXd& v, int rows, int cols, Eigen::MatrixXd& m) {
    m.resize(rows, cols);
    for (int i = 0; i < v.size(); i++) {
        int i_flr = i / cols, i_rmdr = i % cols;
        m(i_flr, i_rmdr) = v(i);
    }
}

void mat_to_vector(const Eigen::MatrixXd& m, Eigen::VectorXd& v)
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> M(m);
    Eigen::Map<Eigen::RowVectorXd> vector_m(M.data(), M.size());
    v = vector_m;

}

// use at your own risk
void print_sparse_mat(const SparseMat& mat, const std::string& name) {
    std::cout << "printing sparse mat : " << name << std::endl;
    std::cout << "row, col: value" << std::endl;

    for (int i = 0; i < mat.outerSize(); i++) {
        for (SparseMat::InnerIterator it(mat, i); it; ++it) {
            std::cout << it.row() << ", " << it.col() << ": " << it.value() << std::endl;
        }

    }
    std::cout << "done with printing sparse mat : " << name << std::endl;
}

int rank_mat(const Eigen::MatrixXd& m)
{
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(m);
    return qr.rank();
}

int rank_mat(const SparseMat& m)
{
    SparseMat m_copy = m;
    m_copy.makeCompressed();
    //SparseMat::makeCompressed(m);
    Eigen::SparseQR<SparseMat, Eigen::COLAMDOrdering<int>> qr(m_copy);
    return qr.rank();
}

void Q_basis_mat(const Eigen::MatrixXd& m, Eigen::MatrixXd& q)
{
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(m);
    q = qr.matrixQ() * Eigen::MatrixXd::Identity(m.rows(), qr.rank());
    return;
}

void Q_basis_mat(const SparseMat& m, SparseMat& q)
{
    //std::cout << "THIS MIGHT BE BUGGED !!" << std::endl;

    SparseMat m_copy = m;
    m_copy.makeCompressed();
    Eigen::SparseQR<SparseMat, Eigen::COLAMDOrdering<int>> qr(m_copy);
    q = qr.matrixQ();
    q.conservativeResize(q.rows(), qr.rank());
    return;
}

double compute_smallest_eigen_value(const SparseMat& m, int n)
{
    using namespace Spectra;
    // Construct matrix operation object using the wrapper class
    SparseSymShiftSolve<double> op(m);

    for (int i = 1; i <= 4; i++) {
        // Construct eigen solver object with shift 0
        // This will find eigenvalues that are closest to 0

        SymEigsShiftSolver<SparseSymShiftSolve<double>> eigs(op, n, n*pow(5,i) , 1e-9);
        //std::cout << " -> factorization with given shift succesful" << std::endl;
        eigs.init();
        eigs.compute(SortRule::LargestMagn);

        if (eigs.info() == CompInfo::Successful)
        {
            Eigen::VectorXd evalues = eigs.eigenvalues();
            //std::cout << "Eigenvalues found: " <<std::endl;
            //for (int j = 0; j < evalues.size(); j++) {
            //    std::cout << j << " " << evalues(j) << std::endl;
            //}
            double aaa = evalues(0);
            for (int j = 1; j < evalues.size(); j++) {
                if (evalues(j) > 0) {
                    if (aaa < 0) {
                        aaa = evalues(j);
                    }
                    else if (evalues(j) < aaa)
                        aaa = evalues(j);
                }
            }
            //return eigs.eigenvalues()(0);
            if (aaa < 0) {
                if (n > 20) return 1;
                return compute_smallest_eigen_value(m, n * 2);
            }
            return aaa;
        }
        /*else {
            std::cout << "Could not find closest eigenvalue closest to 0 ..." << std::endl;
        }*/
    }
    std::cout << "Could not find eigenvalue closest to 0 ..."<<std::endl;
    return 0.0;
}

void convert_int_to_alphas(long long int k, int len_alphas, Eigen::VectorXd& alphas)
{
    alphas = Eigen::VectorXd::Zero(len_alphas);

    int power_idx = len_alphas - 1;

    while (k > 0) {

        if (pow(2, power_idx) <= k) {
            k = k - pow(2, power_idx);
            alphas(len_alphas - 1 - power_idx) = 1;
        }
        power_idx--;
    }
    
    return;
}

void reflection_pt_to_plane(Eigen::Vector3d& p, double A, double B, double C, double D) {
    double d = distance_pt_to_plane(p, A, B, C, B);
    Eigen::Vector3d un;
    normal_to_plane(A, B, C, D, un, true);
    p = p - un * d * 2;
    return;
}

double distance_pt_to_plane(const Eigen::Vector3d& p, double A, double B, double C, double D)
{
    Eigen::Vector3d un, n;
    normal_to_plane(A, B, C, D, un, true);
    normal_to_plane(A, B, C, D, n, false);
    double distance = un.dot(p) + (D / n.norm());
    return distance;
}

void normal_to_plane(double A, double B, double C, double D, Eigen::Vector3d& n, bool unit)
{
    n << A, B, C;
    if (unit) {
        n = n / n.norm();
    }
    return;
}


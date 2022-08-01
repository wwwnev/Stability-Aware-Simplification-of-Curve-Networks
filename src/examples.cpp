#include "examples.h"
#include "surface.h"
#include "load_curves.h"
#include "cn_eigen_solver.h"
#include "cn_greedy_solver.h"
#include "test_network_generator.h"
#include "polyscope/surface_mesh.h"
#include "cn_newton_solver.h"

auto load_results = [](std::string s, int max_i) {

	// loads multiple curve networks and displays them

	for (int i = 0; i <= max_i; i++) {
		std::cout << "i " << i << std::endl;
		curve_network tmp_cn;
		load_curves_for_display(s + std::to_string(i) + ".txt", tmp_cn);
		tmp_cn.update_for_display();
		std::string display_s = "";
		if (max_i < 100) {
			if (i < 10) {
				display_s.append("0");
				
			}
			display_s.append(std::to_string(i));
		}
		else {
			if (i < 10) {
				display_s.append("0");
			}
			if (i < 100) {
				display_s.append("0");
			}
			display_s.append(std::to_string(i));
		}
		tmp_cn.display(display_s, false, false);
	}
	return;
};

void roof(curve_network& cn)
{
	double scale = 70;// 100;
	cn = curve_network();
	std::string fullpath = "input/meshes/roof.obj";
	surface s = surface(fullpath, true);
	s.m_V *= scale;

	// loading the curves from an external file
	load_curves("input/curves/roof_75.txt", cn, scale); // maybe only need load_curves, after some bug fixing
 
	cn.update();
	cn.resample_cn(2 * scale);
	cn.update();
	cn.lock_initial_cn();

	Eigen::VectorXd alphas = Eigen::VectorXd::Ones(cn.m_curves.size());

	bool do_optim_w_relax = true;
	bool do_visualize = false;

	if (do_optim_w_relax) {

		double budget = 0;
		cn.compute_per_curve_length();
		budget = cn.m_curve_len.sum();
		budget = budget * 0.05;

		// here, output files will be located in "output/greedy/relaxed/" and will be named "roof_s70_something.txt"
		std::string name = "roof";

		cn_greedy_solver cn_gs(cn);
		cn_gs.solve_with_relax(alphas, name, budget);

	}
	else if (do_visualize) {

		// if you want to visualize a particular iteration. you can find the "a" vector in the log file -> in "output/greedy/relaxed/XYZ_greedy.txt"

		std::vector<double> a = { 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1 };
		for (int i = 0; i < a.size(); i++)
			a[i] *= 0.15;
		
		// this will display each curve seperately:
		cn.display_per_curve_with_radius("optim ", a);

		// to display as one curve network, build a new one or remove curves from curve_network cn.
	}

	// display the curve network cn in its current state, and the surface.

	cn.display("full");
	polyscope::registerSurfaceMesh("my mesh", s.m_V, s.m_F);
}

void stadium(curve_network& cn, int nb_curves)
{
	std::string fullpath2 = "input/meshes/stadium.obj";
	surface s = surface(fullpath2, true);
	double scale = 30;

	cn = curve_network();
	
	s.m_V *= scale; // scaling the surface, will also scale the curves traced on the surface.

	generate_stadium(cn, s, nb_curves);
	cn.resample_cn(50.0);
	cn.update();
	cn.lock_initial_cn();

	Eigen::VectorXd alphas = Eigen::VectorXd::Ones(cn.m_curves.size());

	bool do_optim_w_relax = true;
	bool do_random_w_relax = false;
	
	// random algo
	if (do_random_w_relax) {
		// relax initial set
		cn_newton_solver n_solver(cn);
		n_solver.solve();
		n_solver.m_CN_current.display("aa");
		cn = n_solver.m_CN_current;

		// get the budget, based on the result from the greedy algorithm
		std::vector<double> a = { 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1 };
		double budget = 0.0;
		cn.compute_per_curve_length();
		for (int i = 0; i < a.size(); i++) {
			if (a[i] > 0)
				budget += cn.m_curve_len[cn.m_curves[i]];
		}
		// time is machine dependant, run this on the same computer
		double time = 1016132;
		
		cn_eigen_solver cn_es(cn, 10);
		Eigen::VectorXd random_alphas;

		cn_es.random_sample_brute_force_v2(1.0, time, budget, random_alphas, {}, "stadium_random", {});
		cn.display_per_curve_with_radius("optim ", random_alphas);
		cn.write_curves("stadium_random", random_alphas);
	}

	else if (do_optim_w_relax) {

		double budget = 0;
		cn.compute_per_curve_length();
		budget = cn.m_curve_len.sum();
		budget = budget * 0.1;

		std::string name = "stadium";
		cn_greedy_solver cn_gs(cn);
		cn_gs.solve_with_relax(alphas, name, budget);
	}

	cn.display("full");
	polyscope::registerSurfaceMesh("my mesh", s.m_V, s.m_F);
}

void shell(curve_network& cn)
{
	double scale = 50.0;
	cn = curve_network();
	std::string fullpath = "input/meshes/shell.obj";
	surface s = surface(fullpath, true);
	s.m_V *= scale;

	load_curves("input/curves/shell_3.txt", cn, s, scale);
	load_curves("input/curves/shell_100.txt", cn, s, scale); // manually removed 40 

	cn.resample_cn(50.0);
	cn.update();
	cn.lock_initial_cn();

	Eigen::VectorXd alphas = Eigen::VectorXd::Ones(cn.m_curves.size());
	
	bool do_optim_w_relax = true;

	if (do_optim_w_relax) {

		double budget = 0;
		cn.compute_per_curve_length();
		budget = cn.m_curve_len.sum();
		budget = budget * 0.05;

		std::string name = "shell";
		cn_greedy_solver cn_gs(cn);
		cn_gs.solve_with_relax(alphas, name, budget);
	}

	cn.display("full");
	polyscope::registerSurfaceMesh("my mesh", s.m_V, s.m_F);
}

void arcshell(curve_network& cn)
{
	double scale = 50.0;
	cn = curve_network();
	std::string fullpath = "input/meshes/shell.obj";
	surface s = surface(fullpath, true);
	s.m_V *= scale;

	load_curves("input/curves/shell_arcs.txt", cn, s);
	load_curves("input/curves/shell_elipse.txt", cn, s);
	load_curves("input/curves/shell_100_for_arcs.txt", cn, s);

	cn.m_curve_set.m_all_curves[6].close_curve();

	// some curves that we want to keep
	std::vector<int> keep;
	for (int i = 0; i < 7; i++) {
		keep.push_back(i);
	}

	cn.resample_cn(2.0);
	cn.update();
	cn.scale(scale);
	cn.update();
	cn.lock_initial_cn();

	Eigen::VectorXd alphas = Eigen::VectorXd::Ones(cn.m_curves.size());

	bool do_optim_w_relax = true;
	
	if (do_optim_w_relax) {

		double budget = 0;
		cn.compute_per_curve_length();
		budget = cn.m_curve_len.sum();
		budget = budget * 0.1;

		std::string name = "arcshell";
		cn_greedy_solver cn_gs(cn);
		cn_gs.solve_with_relax(alphas, name, budget, keep);

	}

	cn.display("full");
	polyscope::registerSurfaceMesh("my mesh", s.m_V, s.m_F);
}

void kagome(curve_network& cn) {
	double scale = 60;
	cn = curve_network();
	std::string fullpath = "input/meshes/kagome.obj";
	surface s = surface(fullpath, true);
	s.m_V *= scale;

	load_curves("input/curves/kagome_114_regular.txt", cn, s);

	cn.resample_cn(2);
	cn.update();
	cn.lock_initial_cn();

	Eigen::VectorXd alphas = Eigen::VectorXd::Ones(cn.m_curves.size());

	bool do_optim_w_relax = true;

	if (do_optim_w_relax) {

		double budget = 0;
		cn.compute_per_curve_length();
		budget = cn.m_curve_len.sum();
		budget = budget * 0.05;

		std::string name = "kagome";
		cn_greedy_solver cn_gs(cn);
		cn_gs.solve_with_relax(alphas, name, budget);
	}

	cn.display("full", true);
	polyscope::registerSurfaceMesh("my mesh", s.m_V, s.m_F);
}

void bunny(curve_network& cn, int nb_curves) {

	double scale = 22.222;
	cn = curve_network();
	std::string fullpath = "input/meshes/bunny_modified.obj";
	surface s = surface(fullpath, true);
	srand(9);
	
	generate_bunny(cn, s, nb_curves);
	s.m_V *= scale;

	cn.resample_cn(2);
	cn.update();
	cn.scale(scale);
	cn.update();
	cn.lock_initial_cn();
	
	Eigen::VectorXd alphas = Eigen::VectorXd::Ones(cn.m_curves.size());

	bool do_optim_w_relax = true;
	if (do_optim_w_relax) {

		double budget = 0;
		cn.compute_per_curve_length();
		budget = cn.m_curve_len.sum();
		budget = budget * 0.1;

		std::string name = "bunny";
		cn_greedy_solver cn_gs(cn);
		std::vector<int> keep = { 60, 65, 67 };
		cn_gs.solve_with_relax(alphas, name, budget, keep);
	}

	cn.display("full");
	polyscope::registerSurfaceMesh("my mesh", s.m_V, s.m_F);
}

void hill(curve_network& cn, int nb_curves)
{
	double scale = 50.0;
	cn = curve_network();
	std::string fullpath = "input/meshes/hill.obj";
	surface s = surface(fullpath, true);
	s.m_V *= scale;

	load_curves("input/curves/hill_3.txt", cn, s, scale);
	generate_hill(cn, s, nb_curves);

	cn.resample_cn(2.0 * scale);
	cn.update();
	cn.lock_initial_cn();

	Eigen::VectorXd alphas = Eigen::VectorXd::Ones(cn.m_curves.size());
	
	bool do_optim_w_relax = true;

	if (do_optim_w_relax) {

		double budget = 0;
		cn.compute_per_curve_length();
		budget = cn.m_curve_len.sum();
		budget = budget * 0.05;

		std::string name = "hill";
		cn_greedy_solver cn_gs(cn);
		std::vector<int> keep = { 0,1 };
		cn_gs.solve_with_relax(alphas, name, budget, keep);

	}

	cn.display("full",false);
	polyscope::registerSurfaceMesh("my mesh", s.m_V, s.m_F);

}

void tent(curve_network& cn, int nb_curves) {

	double scale = 50.0;
	cn = curve_network();
	std::string fullpath = "input/meshes/tent.obj";
	surface s = surface(fullpath, true);

	load_curves("input/curves/tent_3_regular_curves.txt", cn, s);
	load_curves("input/curves/tent_curves.txt", cn, s);

	generate_tent(cn, s, nb_curves);
	s.m_V *= scale;

	cn.resample_cn(2.0);
	cn.update();
	cn.scale(scale);
	cn.update();
	cn.lock_initial_cn();

	Eigen::VectorXd alphas = Eigen::VectorXd::Ones(cn.m_curves.size());

	bool do_optim_w_relax = true;

	if (do_optim_w_relax) {

		double budget = 0;
		cn.compute_per_curve_length();
		budget = cn.m_curve_len.sum();
		budget = budget * 0.05;

		std::string name = "tent";
		cn_greedy_solver cn_gs(cn);
		std::vector<int> keep = { 0,1 };
		cn_gs.solve_with_relax(alphas, name, budget, keep);
	}

	cn.display("full", true);
	polyscope::registerSurfaceMesh("my mesh", s.m_V, s.m_F);
}

void tower(curve_network& cn) {
	double scale = 10.0;
	cn = curve_network();
	std::string fullpath = "input/meshes/tower.obj"; 
	surface s = surface(fullpath, true);
	s.m_V *= scale;

	load_curves("input/curves/tower_top_ring.txt", cn, s);
	load_curves("input/curves/tower_32.txt", cn, s);

	cn.update();
	cn.scale(scale);
	cn.resample_cn(20.0);
	cn.update();
	cn.lock_initial_cn();

	Eigen::VectorXd alphas = Eigen::VectorXd::Ones(cn.m_curves.size());

	bool do_random_w_relax = false;
	bool do_optim_w_relax = true;

	if (do_random_w_relax) {
		// relax initial set
		cn_newton_solver n_solver(cn);
		n_solver.solve();
		n_solver.m_CN_current.display("aa");
		cn = n_solver.m_CN_current;

		// budget and time
		std::vector<double> a = { 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0 };
		double budget = 0.0;
		cn.compute_per_curve_length();
		for (int i = 0; i < a.size(); i++) {
			if (a[i] > 0)
				budget += cn.m_curve_len[cn.m_curves[i]];
		}
		double time = 147791;
		
		cn_eigen_solver cn_es(cn, 10);
		Eigen::VectorXd random_alphas;

		cn_es.random_sample_brute_force_v2(1.0, time, budget, random_alphas, {}, "tower_random", {0}); // keep the first curve (closed ring)
		cn.display_per_curve_with_radius("optim ", random_alphas);
		cn.write_curves("relaxed/tower_random", random_alphas);
	}
	else if (do_optim_w_relax) {

		double budget = 0;
		cn.compute_per_curve_length();
		budget = cn.m_curve_len.sum();
		budget = budget * 0.2;

		std::string name = "tower";
		cn_greedy_solver cn_gs(cn);
		std::vector<int> keep = { 0 };
		cn_gs.solve_with_relax(alphas, name, budget, keep);
	}

	cn.display("full", true);
	polyscope::registerSurfaceMesh("my mesh", s.m_V, s.m_F);
}

void your_example(curve_network& cn, std::string name, std::string in_path, std::string out_path, int nb_curves, double budget)
{
	cn = curve_network();
	surface s = surface(in_path, true);

	double min_x = s.m_V.col(0).minCoeff();
	double max_x = s.m_V.col(0).maxCoeff();
	double range_x = max_x - min_x;

	double min_y = s.m_V.col(1).minCoeff();
	double max_y = s.m_V.col(1).maxCoeff();
	double range_y = max_y - min_y;

	double min_z = s.m_V.col(2).minCoeff();
	double max_z = s.m_V.col(2).maxCoeff();
	double range_z = max_z - min_z;

	double range = range_x > range_y ? range_x : range_y;
	range = range > range_z ? range : range_z;

	double scale = 4000.0 / range;

	s.m_V *= scale;

	generate_your_example(cn, s, nb_curves);
	cn.resample_cn(100);
	cn.update();
	cn.lock_initial_cn();

	Eigen::VectorXd alphas = Eigen::VectorXd::Ones(cn.m_curves.size());

	bool do_optim_w_relax = true;

	if (do_optim_w_relax) {

		double max_budget = 0;
		cn.compute_per_curve_length();
		max_budget = cn.m_curve_len.sum();
		budget = max_budget * budget;

		cn_greedy_solver cn_gs(cn);
		std::vector<int> keep = {};

		cn_gs.solve_with_relax(alphas, name, budget, keep);
	}

	cn.display("full", false);
	polyscope::registerSurfaceMesh("my mesh", s.m_V, s.m_F);
}

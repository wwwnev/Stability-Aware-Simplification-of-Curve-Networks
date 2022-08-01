#include <cstdlib> // for rand()
#include <igl/boundary_loop.h>

#include "test_network_generator.h"


void generate_stadium(curve_network& cn, surface& s, int n)
{
	std::vector<curve> curves;
	curve tmp;
	srand(4);
	n = 100;
	for (int i = 0; i < n; i++) {
		//std::cout << "stadium random geodesic curve " << i << std::endl;
		double rand_t;
		rand_t = (double)i / (double)n + 1 / (2 * (double)n);
		double rand_angle = (double)rand() / RAND_MAX * M_PI * 3 / 5 + M_PI / 5;
		int boundary_choice = i % 2;
		s.generate_geodesic_curve(tmp, rand_t, rand_angle, boundary_choice);
		curves.push_back(tmp);
	}

	for (curve c : curves) {
		cn.add_curve(c);
	}
}

void generate_hill(curve_network& cn, surface& s, int n)
{
	std::vector<curve> curves;
	curve tmp;

	for (int i = 0; i < n; i++) {
		std::cout << "sine disk geodesic curve " << i << std::endl;

		double t = (double)i / (double)n + 1 / (2 * (double)n);

		double rand_angle = (double)rand() / RAND_MAX * M_PI / 2 + M_PI / 4;

		s.generate_geodesic_curve(tmp, t, rand_angle, 0);
		curves.push_back(tmp);

	}

	for (curve c : curves) {
		cn.add_curve(c);
	}
}

void generate_bunny(curve_network& cn, surface& s, int n)
{
	std::vector<curve> curves;
	curve tmp;

	for (int i = 0; i < n; i++) {
		std::cout << "sine disk geodesic curve " << i << std::endl;

		double t = (double)i / (double)n + 1 / (2 * (double)n);

		double rand_angle = (double)rand() / RAND_MAX * M_PI / 2 + M_PI / 4;

		s.generate_geodesic_curve(tmp, t, rand_angle, 0);

		curves.push_back(tmp);
	}

	for (curve c : curves) {
		cn.add_curve(c);
	}
}

void generate_bf_test(curve_network& cn, surface& s, int n, int random_seed, bool do_rand)
{
	srand(random_seed);
	std::vector<curve> curves;
	curve tmp;
	double garbage = rand();

	for (int i = 0; i < n; i++) {
		double t;
		if (do_rand)
			t = (double)rand() / RAND_MAX;
		else
			t = (double)i / (double)n + 1 / (2 * (double)n);

		double rand_angle = (double)rand() / RAND_MAX * M_PI / 2 + M_PI / 4;
		
		s.generate_geodesic_curve(tmp, t, rand_angle, 0);

		curves.push_back(tmp);
	}

	for (curve c : curves) {
		cn.add_curve(c);
	}
}

void generate_tent(curve_network& cn, surface& s, int n)
{
	std::vector<curve> curves;
	curve tmp;
	int skip_count = 0, count = 0, i = 0;
	while(count < n) {
		double t = (double)rand() / RAND_MAX;

		double rand_angle = (double)rand() / RAND_MAX * 2 * M_PI / 3 + M_PI / 6;

		s.generate_geodesic_curve(tmp, t, rand_angle, 0);

		i++;

		bool test1 = abs(tmp.m_V(tmp.m_E(0, 0), 2)) > 5e-3;
		bool test2 = abs(tmp.m_V(tmp.m_E(tmp.m_E.rows() - 1, 1), 2)) > 5e-3;

		if (test1 || test2) {
			skip_count++;
			continue;
		}
		count++;
		curves.push_back(tmp);
	}

	std::vector<std::vector<int>> arcs = {
		{184,182,180,178,176,174,172,170,168,166,164,162,160,158,156,154,152,150,148,146,143,144},
		{930,929,932,934,936,938,940,942,944,946,948,950,952,954,956,958,960,962,964,966,968,970},
		{866,865,868,870,872,874,876,878,880,882,94,92,88,86,84,82,80,78,75,76}
	};

	for (auto a : arcs) {
		tmp = curve();
		for (int i = 0; i < a.size() - 1; i++) {
			tmp.add_edge(s.m_V.row(a[i]).transpose(), s.m_V.row(a[i + 1]).transpose());
		}
	}

	for (curve c : curves) {
		cn.add_curve(c);
	}
}

void generate_your_example(curve_network& cn, surface& s, int n)
{
	std::vector<curve> curves;
	curve tmp;
	int skip_count = 0, count = 0, i = 0;

	int nb_closed_boundaries;
	std::vector<std::vector<int>> L;
	igl::boundary_loop(s.m_F, L);
	nb_closed_boundaries = L.size();

	while (count < n) {
		double t = (double)rand() / RAND_MAX;

		double rand_angle = (double)rand() / RAND_MAX * 2 * M_PI / 3 + M_PI / 6;

		int boundary_start = count % nb_closed_boundaries;

		s.generate_geodesic_curve(tmp, t, rand_angle, 0);

		count++;
		curves.push_back(tmp);
	}

	for (curve c : curves) {
		cn.add_curve(c);
	}
}
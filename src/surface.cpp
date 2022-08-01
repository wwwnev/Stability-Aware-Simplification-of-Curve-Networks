#include "surface.h"

surface::surface(std::string path, bool rewrite = false)
{
	m_path = path;
	igl::readOBJ(path, m_V, m_F);

	if (m_V.rows() == 0 || m_F.rows() == 0)
		std::cout << "WARNING: The mesh you are trying to load has not vertices or no faces. Hint: check path, then file" << std::endl;

	// some meshes imported from Rhino (5, 6 or 7, can't remember) had duplicate vertices... this = fix
	fix_duplicate_vertices();

	igl::edges(m_F, m_E);

	generate_edge_adjacency_matrix();

	if (rewrite) {
		igl::writeOBJ(path, m_V , m_F);
	}

	m_CN = curve_network();
}

void surface::fix_duplicate_vertices() {
	//std::cout << "Fixing duplicate vertices" << std::endl;
	//std::cout << " num vert before : " << m_V.rows() << std::endl;
	for (int i = 0; i < m_V.rows(); i++) {
		for (int j = i + 1; j < m_V.rows(); j++) {
			if (m_V.row(i) == m_V.row(j)) {
				removeRow(m_V, j);
				for (int r = 0; r < m_F.rows(); r++) {
					for (int c = 0; c < m_F.cols(); c++) {
						if (m_F(r, c) == j) m_F(r, c) = i;
						else if (m_F(r, c) > j) m_F(r, c) -= 1; // decrease index by one because we remove a row to the matrix these index refers to.
					}
				}
				j--;
			}
		}
	}
	// std::cout << " num vert after : " << m_V.rows() << std::endl;
	//std::cout << "Fixing duplicate vertices : done" << std::endl;
}

// slow, should try and replace with half edge datastructure
void surface::generate_edge_adjacency_matrix() {
	//std::cout << "generating edge adjacency matrix" << std::endl;
	m_EA.conservativeResize(m_E.rows(), 2);
	for (int i = 0; i < m_E.rows(); i++) {
		m_EA(i, 0) = -1; m_EA(i, 1) = -1;
		int found = 0;
		int e0 = m_E(i, 0), e1 = m_E(i, 1);
		for (int j = 0; j < m_F.rows(); j++) {
			if (found == 2) break;
			for (int k = 0; k < 3; k++) {
				if (m_F(j, k) == e0 && m_F(j, (k + 1) % 3) == e1) {
					m_EA(i, 0) = j;
					found++;
				}
				else if (m_F(j, k) == e1 && m_F(j, (k + 1) % 3) == e0) {
					m_EA(i, 1) = j;
					found++;
				}
			}
		}
		found = 0;

	}
	//std::cout << "generating edge adjacency matrix : done" << std::endl;
}

int surface::edge_index(int start, int end) {
	
	for (int i = 0; i < m_E.rows(); i++) {
		if ((m_E(i, 0) == start && m_E(i, 1) == end) || (m_E(i, 0) == end && m_E(i, 1) == start)) {
			return i;
		}
	}
	return -1;
}


// used to algorithmically trace curves over surfaces
void surface::generate_geodesic_curve(curve& ccc, double t, double angle, int boundary_choice)
{
	//std::cout << "Generating a random geodesic curve" << std::endl;
	curve gc = curve(); // geodesic curve we will build

	// #F by #3 adjacent matrix, the element i,j is the id of the triangle
	//        adjacent to the j edge of triangle i
	// used to navigate from one face to the next. Could be replaced with m_EA that i added afterwards
	Eigen::MatrixXi TT;
	igl::triangle_triangle_adjacency(m_F, TT);

	// get boundary loops 
	std::vector<std::vector<int>> L;

	igl::boundary_loop(m_F, L);

	double boundary_len = 0.0;

	for (int i = 0; i < L[boundary_choice].size() - 1; i++) {
		
		int pa, pb;
		pa = L[boundary_choice][i];
		pb = L[boundary_choice][i + 1];

		int e_idx = edge_index(pa, pb);

		double that_edge_len = (m_V.row(m_E(e_idx, 1)) - m_V.row(m_E(e_idx, 0))).norm();

		boundary_len += that_edge_len;
	}

	bool found_param = false;
	double cur_len = 0;
	double e_len;
	Eigen::VectorXd cur_pos, next_pos;
	if (t == 0 || t == 1) {
		std::cerr << "IF YOU ASK FOR ARC LENGTH PARAM T = 0, or T = 1. THE GEODESIC WILL START ON A VERTEX. -> NOT IMPLEMENTED. MOVING T a bit" << std::endl;
		t == 0 ? t = 0.001 : t = 0.999;
	}
	double target_len = t * boundary_len;
	int boundary_size = L[boundary_choice].size();
	int v1_idx = 0; // -> i0
	int v2_idx = (v1_idx + 1) % boundary_size; // -> i1

	while (!found_param) {
		e_len = (m_V.row(L[boundary_choice][v2_idx]) - m_V.row(L[boundary_choice][v1_idx])).norm();
		if (cur_len + e_len == target_len) {
			std::cout << "we are ending up on a vertex. let's wiggle a little bit." << std::endl;
			target_len -= 0.05 * e_len;
		}
		if (cur_len + e_len > target_len) {
			cur_pos = m_V.row(L[boundary_choice][v1_idx]) + (target_len - cur_len) * (m_V.row(L[boundary_choice][v2_idx]) - m_V.row(L[boundary_choice][v1_idx])) / (m_V.row(L[boundary_choice][v2_idx]) - m_V.row(L[boundary_choice][v1_idx])).norm();
			found_param = true;
		}
		else {
			cur_len += e_len;
			v1_idx = v2_idx;
			v2_idx++;
			v2_idx = v2_idx % boundary_size;
		}
	}

	int p0 = L[boundary_choice][v1_idx];
	int p1 = L[boundary_choice][v2_idx];

	int p2;
	int new_tri = -1, last_tri = -1;

	int edge_idx = edge_index(p0, p1);

	new_tri = m_EA(edge_idx, 0) != -1 ? m_EA(edge_idx, 0) : m_EA(edge_idx, 1);


	// p2 is the vertex of the new tri that is neither p0 or p1 :
	if (m_F(new_tri, 0) != p0 && m_F(new_tri, 0) != p1) {
		p2 = m_F(new_tri, 0);
	}
	else if (m_F(new_tri, 1) != p0 && m_F(new_tri, 1) != p1) {
		p2 = m_F(new_tri, 1);
	}
	else if (m_F(new_tri, 2) != p0 && m_F(new_tri, 2) != p1) {
		p2 = m_F(new_tri, 2);
	}

	// for convenience, we make sure p0 comes before p1 in the triangle vertex list, is it really useful ? Does it ever happen ...
	int idx_p0, idx_p1;
	for (int i = 0; i < 3; i++) {
		if (m_F(new_tri, i) == p0 && m_F(new_tri, mod(i - 1, 3)) == p1) {
			std::cout << "switching p0 and p1" << std::endl;
			int tmp = p1;
			p1 = p0;
			p0 = tmp;
			break;
		}
	}

	// defines the line our geodesic follows on the current plane (il for intersecting line)
	Eigen::Vector3d il0, il1;

	// vector in the direction of the line current intersection point lies on (l10) and in the direction of another line on the current plane
	Eigen::Vector3d l10, l21;

	// normal to the current plane
	Eigen::Vector3d normal;

	// new_point is the correct intersection point, new_point1 and 2 are the 2 possible intersection points, we choose the one that lies within the boundary of the current triangle
	Eigen::Vector3d new_point, new_point1, new_point2;
	bool test1, test2;

	// the pts defining the intersected edge -> they will become p0 and p1 when looping.
	int intersected_p1, intersected_p2;

	bool hit_boundary = false;
	bool start_from_vertex = false; // if we start from a vertex, we are looking for intersection only one edge.

	// first point of our geodesic segment. 
	il0 = cur_pos;

	int count = 0;
	
	while (!hit_boundary) {
		l10 = m_V.row(p1) - m_V.row(p0);
		l21 = m_V.row(p2) - m_V.row(p1);

		normal = l10.cross(l21);
		normal.normalize();
		il1 = il0 + (cos(angle) * l10 + sin(angle) * normal.cross(l10));

		test1 = find_line_intersection(m_V.row(p1), m_V.row(p2), il0, il1, new_point1);

		// we only check with a second edge of the triangle if il0 is not on a vertex of a triangle
		if (!start_from_vertex) {
			test2 = find_line_intersection(m_V.row(p2), m_V.row(p0), il0, il1, new_point2);
		}

		else test2 = false;
		bool tmp_flag = false;
		if (!test1 && !test2) {
			std::cerr << "Cant intersect with either of the 2 remaining triangle edges..." << std::endl;
			return;
		}
		else if (!test1) {
			new_point = new_point2;

			// defining the line that has been intersected, we will rotate this line on the plane of the next triangle.
			intersected_p1 = p2;
			intersected_p2 = p0;
		}
		else if (!test2) {
			new_point = new_point1;
			intersected_p1 = p1;
			intersected_p2 = p2;
		}
		else {
			std::cout << "We intersect both edges, meaning we hit a vertex !" << std::endl;
			std::cout << "intersection pt 1 : " << new_point1.transpose();
			std::cout << " - intersection pt 2 : " << new_point2.transpose() << std::endl << std::endl;
			new_point = m_V.row(p2);

			gc.add_edge(il0, new_point);

			// The case where a vertex is hit is handled inside this else

			// 1. check for boundary hit
			for (auto loop : L)
				for (auto e : loop)
					if (p2 == e) {
						hit_boundary = true;
					}
			if (hit_boundary) {
				//std::cout << "hit a boundary vertex" << std::endl;
				break;
			}

			// 2. Compute sum of angles around vertex
			//std::vector<std::pair<int, double>> tri_and_angle;
			// I want 3 values for each neighbor triangle : face idx, reference vertex, angle
			std::vector<std::tuple<int, int, int, double>> tri_and_angle;
			double sum_angle = 0;
			int np0 = p0, np1; // neighbor to p2 : np0, np1 
			int cur_tri = new_tri;

			// iterate over all 1 ring neighbors of p2 (hit vertex)
			while (true) {
				// np1 is the remaining vertex of current triangle
				for (int i = 0; i < 3; i++) {
					if (m_F(cur_tri, i) != np0 && m_F(cur_tri, i) != p2) {
						np1 = m_F(cur_tri, i);
					}
				}
				Eigen::Vector3d tmp_a = m_V.row(np0) - m_V.row(p2);
				Eigen::Vector3d tmp_b = m_V.row(np1) - m_V.row(p2);

				double tmp_angle = acos(tmp_a.dot(tmp_b) / tmp_a.norm() / tmp_b.norm());
				sum_angle += tmp_angle;
				tri_and_angle.push_back(std::make_tuple(cur_tri, np0, np1, tmp_angle));

				// find next triangle
				int edge_idx = edge_index(p2, np1);
				cur_tri = m_EA(edge_idx, 0) != cur_tri ? m_EA(edge_idx, 0) : m_EA(edge_idx, 1);

				// swap np0, imagine a sweeping motion around vertex p2. np1 will be redefined in the next iteration.
				np0 = np1;

				// if we are back to the first triangle, break
				if (cur_tri == new_tri) break;
			}

			// 3. Find in which triangle is the bisector and with which angle w/r to an edge of that triangle
			double bis_angle = sum_angle / 2.0, sub_angle; // sub_angle : we start measuring angle in the "middle" of the current triangle, at the edge traced by the curve
			sum_angle = 0.0;

			Eigen::Vector3d tmp_a = m_V.row(p0) - m_V.row(p2);
			Eigen::Vector3d tmp_b = il0 - m_V.row(p2).transpose();
			sub_angle = acos(tmp_a.dot(tmp_b) / tmp_a.norm() / tmp_b.norm());

			for (auto e : tri_and_angle) {
				//std::cout<<std::get<0>(e) << " " << std::get<3>(e)<< std::endl;
				sum_angle += std::get<3>(e);
				if (sum_angle == bis_angle) {
					std::cerr << "geodesic curve aligning directly with a mesh edge not yet supported";
				}
				else if (sum_angle - sub_angle > bis_angle) {
					new_tri = std::get<0>(e);
					p0 = p2;
					p1 = std::get<1>(e);
					p2 = std::get<2>(e);
					angle = (bis_angle + sub_angle) - (sum_angle - std::get<3>(e));
					break;
				}
			}
			il0 = new_point;

			// We will jump directly to the next iteration of the 'big' while loop. 
			count++;
			start_from_vertex = true;
			continue;
		}

		count++;
		if (tmp_flag) break;

		// add edge to curve
		gc.add_edge(il0, new_point);

		// new angle. reminder : angle = acos( a . b / |a||b| )
		Eigen::Vector3d a = m_V.row(intersected_p2) - m_V.row(intersected_p1);
		Eigen::Vector3d b = new_point - il0;
		double new_angle = acos(a.dot(b) / a.norm() / b.norm());

		// replace values of angle, last_tri, new_tri, p0, p1, p2, il0
		angle = new_angle;
		p0 = intersected_p1; // intersected line becomes "base" line for next triangle crossed
		p1 = intersected_p2;
		last_tri = new_tri;
		il0 = new_point;

		int edge_idx = edge_index(intersected_p1, intersected_p2);
		new_tri = m_EA(edge_idx, 0) != new_tri ? m_EA(edge_idx, 0) : m_EA(edge_idx, 1);

		// checking for boundary hit
		if (new_tri == -1) {
			//std::cout << "hit boundary !!" << std::endl;
			hit_boundary = true;
			break;
		}

		// p2 is the vertex of the new tri that is neither p0 or p1 :
		if (m_F(new_tri, 0) != p0 && m_F(new_tri, 0) != p1) {
			p2 = m_F(new_tri, 0);
		}
		else if (m_F(new_tri, 1) != p0 && m_F(new_tri, 1) != p1) {
			p2 = m_F(new_tri, 1);
		}
		else if (m_F(new_tri, 2) != p0 && m_F(new_tri, 2) != p1) {
			p2 = m_F(new_tri, 2);
		}
		start_from_vertex = false;
	}

	ccc = gc;
	return;
}

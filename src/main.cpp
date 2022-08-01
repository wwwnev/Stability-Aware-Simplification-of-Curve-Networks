#include <iostream>
#include <string>
#include <sstream>
#include <filesystem>

#include <polyscope/polyscope.h>
#include "polyscope/view.h"

#include "examples.h"

int main(int argc, char* argv[])
{
	// The program can take arguments as input. 

	//std::string fullpath = "input/meshes/1wave.obj"; // relative path, will need to test with .exe .....
	//std::string fullpath = "input/meshes/tube_bridge.obj";
	//int nb_curves = 15;
	//int seed = 2;

	curve_network cn;

	if (argc == 2) {
		int seed = 2;
		srand(seed);

		if (std::string(argv[1]) == "hill") {
			std::cout << "Computing the 'hill' example (fig 1)." << std::endl;
			hill(cn, 70);
		}
		else if (std::string(argv[1]) == "kagome") {
			std::cout << "Computing the 'kagome' example (fig 7a)." << std::endl;
			kagome(cn);
		}
		else if (std::string(argv[1]) == "roof") {
			std::cout << "Computing the 'roof' example (fig 7b)." << std::endl;
			roof(cn);
		}
		else if (std::string(argv[1]) == "stadium") {
			std::cout << "Computing the 'stadium' example (fig 10a)." << std::endl;
			stadium(cn, 100);
		}
		else if (std::string(argv[1]) == "tower") {
			std::cout << "Computing the 'tower' example (fig 10d)." << std::endl;
			tower(cn);
		}
		else if (std::string(argv[1]) == "shell") {
			std::cout << "Computing the 'shell' example (fig 11a)." << std::endl;
			shell(cn);
		}
		else if (std::string(argv[1]) == "bunny") {
			std::cout << "Computing the 'bunny' example (fig 11b)." << std::endl;
			bunny(cn, 75);
		}
		else if (std::string(argv[1]) == "tent") {
			std::cout << "Computing the 'tent' example (fig 11c)." << std::endl;
			tent(cn, 66);
		}
		else if (std::string(argv[1]) == "arcshell") {
			std::cout << "Computing the 'arcced shell' example (fig 11d)." << std::endl;
			arcshell(cn);
		}
		else {
			std::cout << "Please provide a valid example name (hill, kagome, roof, stadium, stadium, tower, shell, bunny, tent, arcshell)." << std::endl;
		}
	}

	else if (argc >= 3)
	{
		std::string name = std::string(argv[1]);
		std::string in_path = "input/meshes/";
		in_path.append(argv[1]);

		std::string out_path = "output/greedy/";
		out_path.append(argv[1]);
		
		int nb_curves = atoi(argv[2]);
		
		double budget;
		if (argc >= 4)
			budget = std::stod(argv[3]);
		else
			budget = 0.5;

		int seed;
		if (argc >= 5)
			seed = atoi(argv[4]);
		else
			seed = 2;

		srand(seed);
		
		std::cout << "Your example : " << in_path << ", nb of curves = " << nb_curves << ", budget = " << budget << ", seed = " << seed << std::endl;

		your_example(cn, name, in_path, out_path, nb_curves, budget);
	}
	else {
		std::cout << "Please provide arguments (see read.me for more information)." << std::endl;
	}
	

	//int beginIdx = fullpath.rfind('/');
	//std::string filename = fullpath.substr(beginIdx + 1);
	//filename.erase(filename.length() - 4);

	//std::string outname = ""; // "output/bruteforce/"
	//outname.append(filename);
	//outname.append("_");
	//outname.append(std::to_string(nb_curves));
	//outname.append("_sample2_seed");
	//outname.append(std::to_string(seed));
	//outname.append("_out.txt");

	polyscope::init();


	// --- Call one of these functions for results ---

	

	//fab(cn_geo, 70); // fig 1
	//kagome(cn_geo); // fig 7a
	//maritime_roof_72(cn_geo); // fig 7b
	//beijing_stadium(cn_geo, 100); // fig 10a
	//tower(cn_geo, 24); // fig 10d
	//boring_bean_gridshell_basic(cn_geo); // fig 11a
	//bunny(cn_geo, 75); // fig 11b
	//boring_bean_gridshell(cn_geo); // fig 11d
	//tent(cn_geo, 66); fig 11c

	// --- SOME VISUALISATION STUFF ---

	polyscope::view::upDir = polyscope::view::UpDir::ZUp;
	polyscope::show();

	return 0;
}
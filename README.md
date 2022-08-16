# Stability-Aware Simplication of Curve Networks

This is a C++ implementation of the greedy algorithm presented in our article [Stability-Aware Simplification of Curve Networks](http://www-labs.iro.umontreal.ca/~bmpix/curve_networks/).

## Prerequisites
This project requires C++17 ([Why?](https://devblogs.microsoft.com/cppblog/using-c17-parallel-algorithms-for-better-performance/)).
We have 4 main dependencies:
- Eigen 3.3.9 (or above)
- Libigl
- Spectra 1.0.0 (or above)
- Polyscope

We assume that Eigen3 is already installed on your system, in a default or preferred location. To get the other 3 dependencies, run the following commands from the root of the repository: 
```sh
mkdir deps/
git submodule update --init --recursive
```
This will download Libigl, Spectra and Polyscope to the `deps/` subdirectory.

## Build, compile & run
Your mileage may vary depending on your compiler, but generally the commands needed to install and run the code in release configuration should look like the following. Start from the root of the project.

```sh
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cd ..
build/bin/Release/curvenet.exe [your arguments here]
```

If you want the release configuration (you probably do) and you are using **Visual Studio C++ compiler**, we need to specify that at build time. The commands will look like following:

```sh
mkdir build
cd build/
cmake ..
cmake --build . --config Release
cd ..
build/bin/Release/curvenet.exe [your arguments here]
```

In the case where Eigen3 cannot be found during the `cmake .` call, you might want to try `cmake . -DEigen3_DIR=$HOME/mypackages/share/eigen3/cmake/` instead, replacing the path with one pointing to the `cmake/` subdirectory in your Eigen3 installation directory.
## Command line arguments

Right now, the compiled project supports two cases with it's argument. First, you can recompute the example shown in the article with 
```sh
build/bin/Release/curvenet.exe keyword
```
with one of the following keywords: hill, roof, kagome, stadium, tower, shell, bunny, tent and arcshell.
The second use case allows you to generate and optimize a curve network on a surface of your choice. Note that the surface must be a `.obj` located in the `input/meshes/` subdirectory.

```sh
build/bin/Release/curvenet.exe your_mesh.obj nb_of_curves [budget] [seed]
```
- `nb_of_curve` is an integer that controls the number of generated curves.
- `budget` is a decimal value between 0 and 1, that controls that total length of curves in the optimized curve network compared to the initial one. By default, it is set to `0.5`.
- `seed` is used to control the random number generator. By default, it is set to `2`.

For example, it might want to try something like this:
```sh
build/bin/Release/curvenet.exe your_mesh.obj 1_wave.obj 20 0.3
```

Logs of the process and results are saved in the `output/greedy/` subdirectory.

## BibTeX

```
@article{Neveu:2022:curvenetworks,
author = {William Neveu and Ivan Puhachov and Bernhard Thomaszewski and Mikhail Bessmeltsev},
title = {Stability-Aware Simplification of Curve Networks},
journal = {ACM Transactions on Graphics},
year = {2022},
doi = {10.1145/3528233.3530711}
} 
```

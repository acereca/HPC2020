#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

#include "../../libs/argparse/argparse.hpp"
#include "../../libs/export/export.hpp"
#include "Grid.hpp"

void relaxation_step(Grid const &old, Grid &next);
void relaxation(Grid &g1, Grid &g2, size_t iterations);

int main(int argc, const char **argv)
{
	// get args
	argparse::ArgumentParser ap;
	ap.addArgument("-s", "--size", 1, false);
	ap.addArgument("-i", "--iterations", 1, false);
	ap.parse(argc, argv);

	size_t size = ap.retrieve<size_t>("size");
	size_t iterations = ap.retrieve<size_t>("iterations");

	// clear output file and write header
	std::ofstream fs;
	fs.open("./out/grid.csv", std::ios_base::trunc | std::ios_base::out);
	fs << "# iteration, x dimension, y dimension, grid data" << std::endl;
	fs.close();

	// stencil can't be inplace -> use two grids and swap between each
	// iteration. Final result will be in g1 if iterations is even and in g2
	// if iterations is odd
	Grid g1(size);
	Grid g2(size);

	auto start = std::chrono::high_resolution_clock::now();
	relaxation(g1, g2, iterations);
	auto stop = std::chrono::high_resolution_clock::now();
	auto time_diff =
	    std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start)
		.count();
	std::cout << size << iterations
		  << time_diff / static_cast<double>(iterations) << std::endl;
#ifndef VERBOSE
	if (iterations % 2) {
		g2.to_csv("./out/grid.csv", iterations - 1);
	} else {
		g1.to_csv("./out/grid.csv", iterations - 1);
	}
#endif /* ifndef  */

	return 0;
}

void relaxation(Grid &g1, Grid &g2, size_t iterations)
{
	for (size_t i = 0; i < iterations; ++i) {
		if (i % 2) {
			relaxation_step(g2, g1);
#ifdef VERBOSE
			g1.to_csv("./out/grid.csv", iterations);
#endif
		} else {
			relaxation_step(g1, g2);
#ifdef VERBOSE
			g2.to_csv("./out/grid.csv", iterations);
#endif
		}
	}
}

void relaxation_step(Grid const &old, Grid &next)
{
	// start loops at 1 and stop one index early since we want to ignore the
	// borders
	constexpr double phi = 24. / 100.;
	for (size_t i = 1; i < old.dimx - 1; ++i) {
		for (size_t j = 1; j < old.dimy - 1; ++j) {
			next(i, j) =
			    old(i, j) + phi * ((-4) * old(i, j) +
			                       old(i + 1, j) + old(i - 1, j) +
			                       old(i, j + 1) + old(i, j - 1));
		}
	}
}

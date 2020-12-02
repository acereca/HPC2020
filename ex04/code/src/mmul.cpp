#include <chrono>
#include <iostream>
#include <random>
#include <string>

#include <argparse.hpp>

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 5120000 / 8
#endif /* ifndef BLOCK_SIZE */

void mmul_unoptimized(double const *A, double const *B, double *C, size_t dim);
void mmul_optimized(double const *A, double const *B, double *C, size_t dim);

int main(int argc, const char **argv) {
	// get args
	argparse::ArgumentParser ap;
	ap.addArgument("-s", "--size", 1, false);
	ap.parse(argc, argv);

	size_t size = ap.retrieve<size_t>("size");
	// setup matrices
	double *A = (double *)calloc(size * size, sizeof(double));
	double *B = (double *)calloc(size * size, sizeof(double));
	double *C = (double *)calloc(size * size, sizeof(double));
	// setup rng
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(0.0, 100.0);
	// fill with random values
	for (size_t i = 0; i < size * size; ++i) {
		A[i] = dist(gen);
	}
	for (size_t i = 0; i < size * size; ++i) {
		B[i] = dist(gen);
	}

	auto start_unopt = std::chrono::high_resolution_clock::now();
	mmul_unoptimized(A, B, C, size);
	auto stop_unopt = std::chrono::high_resolution_clock::now();
	auto time_unopt = std::chrono::duration_cast<std::chrono::nanoseconds>(
			      stop_unopt - start_unopt)
			      .count();

	auto start_opt = std::chrono::high_resolution_clock::now();
	mmul_optimized(A, B, C, size);
	auto stop_opt = std::chrono::high_resolution_clock::now();
	auto time_opt = std::chrono::duration_cast<std::chrono::nanoseconds>(
			      stop_opt - start_opt)
			      .count();

	double FLOP = size * size * size * 2.0;
	std::cout << "======================" << std::endl;
	std::cout << "Unoptimized :" << std::endl;
	std::cout << "Time " << time_unopt << " ns" << std::endl;
	std::cout << "=> " << (FLOP / time_unopt) << "GFLOP/s"
		  << std::endl;
	std::cout << "----------------------" << std::endl;
	std::cout << "Optimized :" << std::endl;
	std::cout << "Time " << time_opt << " ns" << std::endl;
	std::cout << "=> " << (FLOP / time_opt) << "GFLOP/s"
		  << std::endl;
	std::cout << "======================" << std::endl;

	free(A);
	free(B);
	free(C);

	return 0;
}

void mmul_unoptimized(double const *A, double const *B, double *C, size_t dim) {
	/*
	 * Naive matrix multiply for square matrices (NxN)
	 * A, B: input matrices
	 * C: output matrix (has to be initialized to zero)
	 * dim: number of rows (or columns in this case) of the matrix
	 */
	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = 0; j < dim; ++j) {
			for (size_t k = 0; k < dim; ++k) {
				C[i * dim + j] +=
				    A[i * dim + k] * B[k * dim + j];
			}
		}
	}
}

void mmul_optimized(double const *A, double const *B, double *C, size_t dim) {
	/*
	 * Optimized version of matrix multiply for square matrices (NxN)
	 * Memory bandwith is increased by better cache utilization via
	 * block wise access to elements instead of linear traversal
	 * A, B : input matrices
	 * C: output matrix (has to be initialized to zero)
	 * dim: number of rows (or columns in this case) of the matrix
	 */
    #pragma omp parallel for collapse(3) shared(C)
	for (size_t ii = 0; ii < dim; ii += BLOCK_SIZE) {
		for (size_t jj = 0; jj < dim; jj += BLOCK_SIZE) {
			for (size_t kk = 0; kk < dim; kk += BLOCK_SIZE) {
				for (size_t i = ii; i < std::min(ii + BLOCK_SIZE, dim); ++i) {
					for (size_t j = jj; j < std::min(jj + BLOCK_SIZE, dim); ++j) {
						for (size_t k = kk; k < std::min(kk + BLOCK_SIZE, dim); ++k) {
							C[i * dim + j] +=
							    A[i * dim + k] * B[k * dim + j];
						}
					}
				}
			}
		}
	}
}

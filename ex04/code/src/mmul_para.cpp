#include <mpi.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <argparse.hpp>

#include "debug.hpp"

void multiply_parallel(double const *A, double const *B, double *c, size_t mdim);
void multiply_parallel(std::vector<double> &A, std::vector<double> &B, std::vector<double> &C, size_t mdim);
void transpose(std::vector<double> &M);
void print_matrix(double const *M, size_t mdim);
void print_matrix(std::vector<double> &M, size_t mdim);

int main(int argc, char const *argv[]) {

    // Command line arguments
    argparse::ArgumentParser ap;
    ap.addArgument("-n", "--matrix-dim", 1, false);
    ap.addArgument("-p", "--printing", 1);
    ap.parse(argc, argv);

    // MPI_Init(&argc, &argv);
    MPI::Init();

    size_t mdim = ap.retrieve<size_t>("matrix-dim");
    int processes, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_DBG::enable_debugging(rank);

    size_t blocksize = ((mdim+processes-1)/processes)*mdim;

    // row-column order
    // double *A = (double *) calloc(mdim * mdim, sizeof(double));
    std::vector<double> A;
    A.resize(mdim*mdim, 0);
    // double *B = (double *) calloc(mdim * mdim, sizeof(double));
    std::vector<double> B;
    B.resize(mdim*mdim, 0);

    double *C = (double *) calloc(blocksize*processes, sizeof(double));
    // double *cinterm = (double *) calloc(
        // blocksize,
        // sizeof(double)
    // );
    std::vector<double> cinterm;
    cinterm.resize(blocksize, 0);

    for (size_t i = 0; i < mdim; i++) {
        for (size_t j = 0; j < mdim; j++) {
            A[i*mdim+j] = i+j;
            B[i*mdim+j] = i*j;
        }
    }

    transpose(B);


    auto starttime = std::chrono::high_resolution_clock::now();
    multiply_parallel(A, B, cinterm, mdim);

    // std::cout << "rank " << rank << ": ";
    // for (size_t i = 0; i < blocksize; i++)
        // std::cout << cinterm[i] << ", ";
    // std::cout << std::endl;

    MPI_Gather(&cinterm[0], blocksize, MPI_DOUBLE, C, blocksize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    auto stoptime = std::chrono::high_resolution_clock::now();
    // auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(stoptime - starttime).count();
    // double flops = mdim * mdim * mdim * 2;

    if (rank == 0){
        // std::cout << "# t/ms, dim, nodes, GFLOPS/s\n";
        // std::cout
            // << diff << ", "
            // << mdim << ", "
            // << processes << ", "
            // << flops/diff/1e6 << "\n";
    }
    if (rank == 0 && ap.retrieve<bool>("printing")) {
        print_matrix(A, mdim);
        print_matrix(B, mdim);
        print_matrix(C, mdim);
    }

    free(C);

    MPI::Finalize();
    return 0;
}

/* worker function
*
*  @param B is transposed
*/
void multiply_parallel(double const *A, double const *B, double *C, size_t mdim) {
    int processes, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ceiling division
    size_t blocksize = ((mdim+processes-1)/processes);
    size_t offset_A = rank*blocksize;

    for (size_t i = 0; i < blocksize; i++) {
        for (size_t j = 0; j < mdim; j++) {
            for (size_t s = 0; s < mdim; s++) {
                C[i*mdim+j] += A[(offset_A+i*mdim)+s] * B[(j*mdim)+s];
            }
        }
    }

}

void multiply_parallel(std::vector<double> &A, std::vector<double> &B, std::vector<double> &C, size_t mdim) {
    int processes, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ceiling division
    size_t blocksize = ((mdim+processes-1)/processes);
    size_t offset_A = rank*blocksize;

    for (size_t i = 0; i < blocksize; i++) {
        for (size_t j = 0; j < mdim; j++) {
            for (size_t s = 0; s < mdim; s++) {
                C.at(i*mdim+j) += A.at((offset_A+i*mdim)+s) * B.at((j*mdim)+s);
            }
        }
    }

}
void transpose(std::vector<double> &M) {
    size_t mdim = M.size();
    for (size_t i = 0; i < mdim; i++) {
        for (size_t j = 0; j < mdim; j++) {
            size_t row = i + mdim * j;
            size_t col = j + mdim * i;
            double tmp = M[row];
            M[row] = M[col];
            M[col] = tmp;
        }
    }
}

void print_matrix(double const *M, size_t mdim) {
    std::cout << "Matrix:" << std::endl;
    for (size_t i = 0; i < mdim; ++i) {
        for (size_t j = 0; j < mdim; ++j) {
            std::cout << M[i * mdim + j] << ", ";
        }
        std::cout << std::endl;
    }
}

void print_matrix(std::vector<double> &M, size_t mdim) {
    std::cout << "# Matrix:" << std::endl;
    for (size_t i = 0; i < mdim; ++i) {
        for (size_t j = 0; j < mdim; ++j) {
            std::cout << M[i * mdim + j] << ", ";
        }
        std::cout << std::endl;
    }
}

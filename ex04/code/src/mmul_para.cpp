#include <mpi.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <argparse.hpp>

#include "debug.hpp"

void multiply_parallel(std::vector<double> const &A, std::vector<double> const &B, std::vector<double> &C, size_t mdim);
void transpose(std::vector<double> &M, size_t mdim);
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

    // size of blocks of A
    size_t blocksize = ((mdim+processes-1)/processes)*mdim;

    // row-column order
    std::vector<double> A;
    A.resize(mdim*mdim, 0);
    std::vector<double> B;
    B.resize(mdim*mdim, 0);

    // double *C = (double *) calloc(blocksize*processes, sizeof(double));
    std::vector<double> C;
    C.resize(blocksize * processes, 0);
    std::vector<double> cinterm;
    cinterm.resize(blocksize, 0);

    for (size_t i = 0; i < mdim; i++) {
        for (size_t j = 0; j < mdim; j++) {
            A[i*mdim+j] = i+j;
            B[i*mdim+j] = i*j;
        }
    }

    transpose(B, mdim);

    auto starttime = std::chrono::high_resolution_clock::now();
    multiply_parallel(A, B, cinterm, mdim);

    MPI_Gather(
        &cinterm[0], blocksize, MPI_DOUBLE,
        &C[0], blocksize, MPI_DOUBLE,
        0, MPI_COMM_WORLD);

    auto stoptime = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(stoptime - starttime).count();
    double flops = mdim * mdim * mdim * 2;

    if (rank == 0){
        std::cout << "# t/ns, dim, nodes, GFLOPS64/s\n";
        std::cout
            << diff << ", "
            << mdim << ", "
            << processes << ", "
            << flops/diff << "\n";
    }
    if (rank == 0 && ap.retrieve<bool>("printing")) {
        print_matrix(A, mdim);
        print_matrix(B, mdim);
        print_matrix(C, mdim);
    }

    MPI::Finalize();
    return 0;
}

/* worker function
*  multiplies part of A and B, writes to C
*
*  @param[in] A - mdim x mdim Matrix
*  @param[in] B - mdim x mdim Matrix (transposed)
*  @param[out] C - blocksize x mdim Matrix output block
*/
void multiply_parallel(std::vector<double> const &A, std::vector<double> const &B, std::vector<double> &C, size_t mdim) {
    int processes, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ceiling division
    size_t blocksize = ((mdim+processes-1)/processes);
    size_t offset_A = rank*blocksize;

    for (size_t i = 0; i < blocksize; i++) {
        for (size_t j = 0; j < mdim; j++) {
            for (size_t s = 0; s < mdim; s++) {
                C[(i*mdim+j)] += A[(offset_A+i*mdim)+s] * B[(j*mdim)+s];
            }
        }
    }

}

void transpose(std::vector<double> &M, size_t mdim) {
    for (size_t i = 0; i < mdim; i++) {
        for (size_t j = i+1; j < mdim; j++) {
            size_t row = i + mdim * j;
            size_t col = j + mdim * i;
            double tmp = M[row];
            M[row] = M[col];
            M[col] = tmp;
        }
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

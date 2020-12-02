#include <mpi.h>
#include <iostream>
#include <argparse.hpp>

#include <debug.hpp>

void multiply_parallel(double const *A, double const *B, double *c, size_t mdim);
void transpose(double *M, size_t mdim);
void print_matrix(double const *M, size_t mdim);

int main(int argc, char const *argv[]) {

    MPI_Init(NULL, NULL);
    // Command line arguments
    argparse::ArgumentParser ap;
    ap.addArgument("-n", "--matrix-dim", 1, false);
    ap.parse(argc, argv);

    size_t mdim = ap.retrieve<size_t>("matrix-dim");
    int processes, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_DBG::enable_debugging(rank);

    size_t blocksize = ((mdim*mdim+processes-1)/processes);

    // row-column order
    double *A = (double *) calloc(mdim * mdim, sizeof(double));
    double *B = (double *) calloc(mdim * mdim, sizeof(double));
    double *C = (double *) calloc(blocksize*processes, sizeof(double));
    double *cinterm = (double *) calloc(
        blocksize,
        sizeof(double)
    );

    for (size_t i = 0; i < mdim; i++) {
        for (size_t j = 0; j < mdim; j++) {
            A[i*mdim+j] = i+j;
            B[i*mdim+j] = i*j;
        }
    }
    transpose(B, mdim);

    multiply_parallel(A, B, cinterm, mdim);

    // std::cout << "rank " << rank << ": ";
    // for (size_t i = 0; i < blocksize; i++)
        // std::cout << cinterm[i] << ", ";
    // std::cout << std::endl;

    // MPI_Gather(cinterm, blocksize, MPI_DOUBLE, C, blocksize, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if (rank == 0) {
        print_matrix(A, mdim);
        print_matrix(B, mdim);
        print_matrix(C, mdim);
    }

    free(A);
    free(B);
    free(C);
    free(cinterm);

    MPI_Finalize();
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

void transpose(double *M, size_t mdim) {
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

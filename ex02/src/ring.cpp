#include <string>
#include <iostream>
#include <chrono>

#include <mpi.h>
#include <argparse.hpp>

int main(int argc,  const char** argv) {

    // Command line arguments
    argparse::ArgumentParser ap;
    ap.addArgument("-n", "--nummsg", 1, false);
    ap.parse(argc, argv);

    // Setup
    MPI::Init();

    int comm_size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) std::cout << "rank;nproc;nummsg;t_total/ms;t_msg/ns\n";

    // Compute
    uint8_t data = 0;
    MPI_Status status;
    int nummsg = ap.retrieve<int>("nummsg");
    auto starttime = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < nummsg; i++) {
        MPI_Send(&data, 0, MPI_BYTE, (rank+1)%comm_size, 0, MPI_COMM_WORLD);
        MPI_Recv(&data, 1, MPI_BYTE, (rank-1)%comm_size, 0, MPI_COMM_WORLD, &status);
    }
    auto endtime = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(endtime-starttime).count();
    std::cout
        << rank << ";"
        << comm_size << ";"
        << nummsg << ";"
        << diff/1000000 << ";"
        << diff/nummsg << "\n";

    // Finalize
    MPI::Finalize();

    return 0;
}

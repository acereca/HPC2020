#include <string>
#include <iostream>
#include <chrono>

#include <mpi.h>
#include <argparse.hpp>

int main(int argc,  const char** argv) {

    // Command line arguments
    argparse::ArgumentParser ap;
    ap.addArgument("-n", "--nummsg", 1, false);
    ap.addArgument("-s", "--msgsize", 1, false);
    ap.addArgument("-N", "--nodes", 1, false);
    ap.parse(argc, argv);

    // Setup
    MPI::Init();

    int comm_size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    // Compute
    MPI_Status status;
    int nummsg = ap.retrieve<int>("nummsg");
    int msgsize = ap.retrieve<int>("msgsize");
    int nodes = ap.retrieve<int>("nodes");
    std::vector<uint8_t> data {};
    data.resize(msgsize*1024, 0);

    if (rank == 0){
        std::cout << "# t_msg/ns, s_msg/kB, nodes\n";
        for (int i = 0; i < nummsg; i++) {
            auto starttime = std::chrono::high_resolution_clock::now();
            MPI_Send(&data[0], msgsize*1024, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&data[0], msgsize*1024, MPI_BYTE, 1, 0, MPI_COMM_WORLD, &status);
            auto endtime = std::chrono::high_resolution_clock::now();
            auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(endtime-starttime).count();
            std::cout
                << diff << ", "
                << msgsize << ", "
                << nodes << "\n";


        }
    } else {
        for (int i = 0; i < nummsg; i++) {
            MPI_Recv(&data[0], msgsize*1024, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Send(&data[0], msgsize*1024, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
        }
    }

    // Finalize
    MPI::Finalize();

    return 0;
}

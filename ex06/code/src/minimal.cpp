#include <mpi.h>
#include <random>
#include <string>
#include <vector>

#include "../../libs/argparse/argparse.hpp"
#include "../../libs/export/export.hpp"
#include "Grid.hpp"
int main(int argc, const char **argv) {
  // get args
  argparse::ArgumentParser ap;
  ap.addArgument("-s", "--size", 1, false);
  ap.addArgument("-i", "--iterations", 1, false);
  ap.parse(argc, argv);

  size_t size = ap.retrieve<size_t>("size");
  size_t iterations = ap.retrieve<size_t>("iterations");

  int processes, rank;

  MPI::Init();
  MPI_Comm WORLD_1D;
  int dims[] = {processes};
  int periods[] = {0};

  MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 0, &WORLD_1D);
  MPI_Finalize();
}

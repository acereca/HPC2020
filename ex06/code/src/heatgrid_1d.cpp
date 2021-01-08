#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <mpi.h>
#include <random>
#include <string>
#include <vector>

#include "../../libs/argparse/argparse.hpp"
#include "../../libs/export/export.hpp"
#include "../../libs/mpi_debug/debug.hpp"
#include "Grid.hpp"

void gather_and_print(size_t size, size_t processes, std::vector<Grid> &g,
                      size_t iteration, MPI_Comm comm);
void relaxation_step(Grid const &old, Grid &next);
// void relaxation(Grid &g1, Grid &g2, size_t iterations);

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
  MPI_Comm_size(MPI_COMM_WORLD, &processes);

  MPI_Comm WORLD_1D;
  const int dims[] = {processes};
  const int periods[] = {false};

  MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, false, &WORLD_1D);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_DBG::enable_debugging(rank);

  // 64 rows / 3 processes -> 22+2 (x3)
  // row 0 and row 23 contain halo elements
  Grid g1(size, (size + processes - 1) / processes + 2);
  Grid g2(size, (size + processes - 1) / processes + 2);
  std::vector<Grid> g = {g1, g2};

  int coord[1];
  MPI_Cart_coords(WORLD_1D, rank, 1, coord);

  if (coord[0] == 0) { // fill top with start condition
    size_t low = g[0].dimx / 4;
    size_t high = 3 * g[0].dimx / 4;
    for (size_t i = low; i < high; ++i) {
      g[0](i, 0) = 127.;
      g[1](i, 0) = 127.;
    }
  }

  std::vector<int> coords = {coord[0] - 1, coord[0] + 1};
  // std::vector<int> neighbors(2, MPI_PROC_NULL);
  int north = MPI_PROC_NULL;
  int south = MPI_PROC_NULL;
  MPI_Cart_shift(WORLD_1D, 0, -1, &rank, &north);
  MPI_Cart_shift(WORLD_1D, 0, 1, &rank, &south);

  auto start = std::chrono::high_resolution_clock::now();
  for (size_t step = 0; step < iterations; step++) {
    // std::cout << coord[0] << ":" << step << std::endl;
    size_t old = step % 2;
    size_t next = (step + 1) % 2;
    relaxation_step(g[old], g[next]);

    // Exchange
    std::cout << "hi from " << coord[0] << "(" << step << ")" << std::endl;
    MPI_Barrier(WORLD_1D);
    MPI_Sendrecv(&g[next](0, 1), size, MPI_DOUBLE, north, 0, &g[next](0, 0),
                 size, MPI_DOUBLE, north, 0, WORLD_1D, NULL);
    MPI_Sendrecv(&g[next](0, g[next].dimy - 2), size, MPI_DOUBLE, south, 0,
                 &g[next](0, g[next].dimy - 1), size, MPI_DOUBLE, south, 0,
                 WORLD_1D, NULL);
    // MPI_Send(&g[next](0, 1), size, MPI_DOUBLE, north, 0, WORLD_1D);
    // MPI_Send(&g[next](0, g[next].dimy - 2), size, MPI_DOUBLE, south, 0,
    // WORLD_1D);
    // MPI_Recv(&g[next](0, 0), size, MPI_DOUBLE, north, 0, WORLD_1D, NULL);
    // MPI_Recv(&g[next](0, g[next].dimy - 1), size, MPI_DOUBLE, south, 0,
    // WORLD_1D, NULL);
    std::cout << "bye from " << coord[0] << "(" << step << ")" << std::endl;
#ifdef VERBOSE
    gather_and_print(size, processes, g, step, MPI_COMM_WORLD);
    // g1.to_csv("./out/grid.csv", iterations);
#endif
  }
  auto stop = std::chrono::high_resolution_clock::now();

  if (coord[0] == 0) {
    std::ofstream fs;
    fs.open("./out/grid.csv", std::ios::out | std::ios::app);
    auto time_diff =
        std::chrono::duration_cast<std::chrono::microseconds>(stop - start)
            .count();
    fs << size << ";" << processes << ";" << iterations << ";" << time_diff
       << "\n";

    fs.close();
  }
  // recvs both successfull

  MPI_Finalize();
}

void gather_and_print(size_t size, size_t processes, std::vector<Grid> &g,
                      size_t iteration, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Grid gall(size, 1 + ((size + processes - 1) / processes) * processes);
  size_t low = gall.dimx / 4;
  size_t high = 3 * gall.dimx / 4;
  for (size_t i = low; i < high; ++i) {
    gall(i, 0) = 127.;
    gall(i, 0) = 127.;
  }
  MPI_Gather(&g[(iteration) % 2](0, 1),
             g[iteration % 2].dimx * (g[iteration % 2].dimy - 2), MPI_DOUBLE,
             &gall(0, 1), g[iteration % 2].dimx * (g[iteration % 2].dimy - 2),
             MPI_DOUBLE, 0, comm);
  if (rank == 0) {
    gall.dimy = size; // to print only 'necessary' data
    gall.to_csv("./out/gall.csv", iteration);
  }
}

void relaxation(Grid &g1, Grid &g2, size_t iterations) {
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

void relaxation_step(Grid const &old, Grid &next) {
  // start loops at 1 and stop one index early since we want to ignore the
  // borders
  constexpr double phi = 24. / 100.;
  for (size_t i = 1; i < old.dimx - 1; ++i) {
    for (size_t j = 1; j < old.dimy - 1;
         ++j) { // dont iterate througb halo elements
      next(i, j) =
          old(i, j) + phi * ((-4) * old(i, j) + old(i + 1, j) + old(i - 1, j) +
                             old(i, j + 1) + old(i, j - 1));
    }
  }
}

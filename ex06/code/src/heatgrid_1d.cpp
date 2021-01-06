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
  std::vector<int> neighbors(2, -1);
  for (size_t i = 0; i < coords.size(); i++) {
    if (coord[i] >= 0 && coord[i] < processes) {
      int c[] = {coords[0]};
      // MPI_Cart_rank(WORLD_1D, c, &neighbors[i]);
      MPI_Cart_shift(WORLD_1D, 0, -1 + 2 * i, &rank, &neighbors[i]);
    }
  }

  auto start = std::chrono::high_resolution_clock::now();
  for (size_t step = 0; step < iterations; step++) {
    std::cout << step << std::endl;
    size_t old = step % 2;
    size_t next = (step + 1) % 2;
    relaxation_step(g[old], g[next]);

    // Exchange
    if (neighbors[0] != -1 || neighbors[0] != MPI_PROC_NULL) {
      MPI_Send(&g[next](0, 1), size, MPI_DOUBLE, neighbors[0], 0, WORLD_1D);
      MPI_Recv(&g[next](0, 0), size, MPI_DOUBLE, neighbors[0], 0, WORLD_1D,
               NULL);
    }

    if (neighbors[1] != -1 || neighbors[1] != MPI_PROC_NULL) {
      MPI_Send(&g[next](0, g[next].dimy - 2), size, MPI_DOUBLE, neighbors[1], 0,
               WORLD_1D);

      MPI_Recv(&g[next](0, g[next].dimy - 1), size, MPI_DOUBLE, neighbors[1], 0,
               WORLD_1D, NULL);
    }
  }
  auto stop = std::chrono::high_resolution_clock::now();
  std::ofstream fs;
  fs.open("./out/grid.csv");
  auto time_diff =
      std::chrono::duration_cast<std::chrono::milliseconds>(stop - start)
          .count();
  fs << size << ";" << processes << ";" << iterations << ";" << time_diff;

  fs.close();
  // recvs both successfull
  MPI_Finalize();
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

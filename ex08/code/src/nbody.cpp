#include <mpi.h>
#include <stdlib.h>

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "../../../libs/argparse/argparse.hpp"
#include "../../../libs/mpi_debug/debug.hpp"

int DELTA_T = 100;		     // s
constexpr double GAMMA = 6.673e-11;  // N m^2 kg^-2
constexpr double EPSILON = 1e-50;    // softening factor
// Will be used to obtain a seed for the random number engine
std::random_device rd;
// Standard mersenne_twister_engine seeded with rd()
std::mt19937 gen(rd());
std::uniform_real_distribution<> mass_distrib(1e24, 9e24);
std::uniform_real_distribution<> space_distrib(0, 1e11);

struct bbuffer {
	std::vector<double> data;

	double *mass;
	double *posX;
	double *posY;
	double *posZ;

	bbuffer() = delete;
	bbuffer(size_t n) {
		data.resize(4 * n, 0);
		mass = &data[0];
		posX = &data[1 * n];
		posY = &data[2 * n];
		posZ = &data[3 * n];
	}
};

struct bodies {
	// contains all data for all bodies of this rank
	// [mass_0..mass_N, posX_0..posX_N, ]
	// double *data;
	std::vector<double> data;

	// masses for all bodies in this rank
	double *mass;

	// same for positions
	double *posX;
	double *posY;
	double *posZ;

	// and velocities
	double *velX;
	double *velY;
	double *velZ;

	size_t n;

	bodies() = delete;
	bodies(size_t n, size_t nulled = 0) : n(n) {
		// data = static_cast<double *>(calloc(7 * n, sizeof(double)));
		data.resize(7 * n, 0);
		mass = &data[0];

		// pointers to sub (we now this is hacky,
		// would require views from c++17, we are stuck with c++11)
		posX = &data[1 * n];
		posY = &data[2 * n];
		posZ = &data[3 * n];
		velX = &data[4 * n];
		velY = &data[5 * n];
		velZ = &data[6 * n];

		for (size_t i = 0; i < (n - nulled); i++) {
			mass[i] = mass_distrib(gen);
			posX[i] = space_distrib(gen);
			posY[i] = space_distrib(gen);
			posZ[i] = space_distrib(gen);
		}
	}
	void print_to_file(std::string fd, size_t step) {
		std::fstream fs;
		if (step == 0) {
			fs = std::fstream(fd, fs.out | fs.trunc);
		} else {
			fs = std::fstream(fd, fs.out | fs.app);
		}
		if (!fs.is_open())
			std::fprintf(stderr, "failed to open file %s\n",
				     fd.c_str());
		else {
			fs << "#step; "
			   << "mass [kg]; "
			   << "pos.x; "
			   << "pos.y; "
			   << "pos.z;" << std::endl;
			for (size_t i = 0; i < n; ++i) {
				fs << step << ", " << mass[i] << ", " << posX[i]
				   << ", " << posY[i] << ", " << posZ[0]
				   << "\n";
			}
		}
		fs.close();
	}
};

size_t ceil(size_t a, size_t b) { return (a + b - 1) / b; }

struct leap_frog_iterator {
	// size_t step_ring = 0;
	const size_t size_ring;
	const size_t n;
	bodies const &b;

	std::vector<std::array<double, 3>> acc;

	leap_frog_iterator() = delete;
	leap_frog_iterator(size_t size_ring, size_t n, bodies &b)
	    : size_ring(size_ring), n(n), b(b) {
		acc.resize(n, {0, 0, 0});
	};

	void update() {
		for (size_t i = 0; i < n; i++) {
			b.posX[i] += b.velX[i] * DELTA_T;
			b.posY[i] += b.velY[i] * DELTA_T;
			b.posZ[i] += b.velZ[i] * DELTA_T;

			b.velX[i] += acc[i][0] * DELTA_T;
			b.velY[i] += acc[i][1] * DELTA_T;
			b.velZ[i] += acc[i][2] * DELTA_T;
		}
		acc.clear();
		acc.resize(n, {0, 0, 0});
	}

	void operator()(bbuffer &msg_objs) {
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				double rX = (b.posX[i] - msg_objs.posX[j]);
				double rY = (b.posY[i] - msg_objs.posY[j]);
				double rZ = (b.posZ[i] - msg_objs.posZ[j]);
				double dist_sq =
				    rX * rX + rY * rY + rZ * rZ + EPSILON;
				double distX = rX / std::sqrt(dist_sq);
				double distY = rY / std::sqrt(dist_sq);
				double distZ = rZ / std::sqrt(dist_sq);
				acc[i][0] =
				    -(GAMMA * msg_objs.mass[j] * distX) /
				    dist_sq;
				acc[i][1] =
				    -(GAMMA * msg_objs.mass[j] * distY) /
				    dist_sq;
				acc[i][2] =
				    -(GAMMA * msg_objs.mass[j] * distZ) /
				    dist_sq;
			}
		}
	}
};

int main(int argc, char const *argv[]) {
	// init stuff
	argparse::ArgumentParser ap;
	ap.addArgument("-n", "--num-objects", 1, false);
	ap.addArgument("-i", "--iterations", 1, false);
	ap.addArgument("-t", "--tstep", 1, false);
	ap.parse(argc, argv);

	size_t size = ap.retrieve<size_t>("num-objects");
	size_t iterations = ap.retrieve<size_t>("iterations");
	DELTA_T = ap.retrieve<size_t>("tstep");

	// int processes = 1;
	// int rank = 0;
	// int prev = 0;
	// int next = 0;
	int processes;
	int rank;
	int prev;
	int next;

	MPI::Init();
	MPI_Comm_size(MPI_COMM_WORLD, &processes);

	MPI_Comm WORLD_1D;
	const int dims[] = {processes};
	const int periods[] = {true};

	MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, false, &WORLD_1D);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		MPI_DBG::enable_debugging(rank);
	}
	int coord[1];
	MPI_Cart_coords(WORLD_1D, rank, 1, coord);

	size_t bodies_per_rank = ceil(size, processes);

	bodies objects(bodies_per_rank);
	std::array<bbuffer, 2> bbuffs{bodies_per_rank, bodies_per_rank};

	if (coord[0] == (processes - 1)) {
		objects = bodies(
		    bodies_per_rank,
		    (bodies_per_rank - (size % processes)) % bodies_per_rank);
	}

	leap_frog_iterator lfi(processes, bodies_per_rank, objects);

	auto t_start = std::chrono::high_resolution_clock::now();
	// compute section
	// MPI_irecv ->0

	MPI_Cart_shift(WORLD_1D, 0, -1, &rank, &prev);
	MPI_Cart_shift(WORLD_1D, 0, 1, &rank, &next);

	MPI_Request sig_recv, sig_send;

	MPI_Irecv(&bbuffs[0].data[0], bodies_per_rank * 4, MPI_DOUBLE, prev, 0,
		  WORLD_1D, &sig_recv);
	MPI_Isend(&objects.data[0], bodies_per_rank * 4, MPI_DOUBLE, next, 0,
		  WORLD_1D, &sig_send);

	for (size_t time_step = 0; time_step < iterations; time_step++) {
		// for processes-1
		for (size_t loop_step = 0;
		     loop_step < static_cast<size_t>(processes - 1);
		     loop_step++) {
			// MPI_wait
			MPI_Wait(&sig_recv, NULL);
			// MPI_isend 0->
			MPI_Wait(&sig_send, NULL);
			MPI_Isend(&bbuffs[loop_step % 2].data[0],
				  bodies_per_rank * 4, MPI_DOUBLE, next, 0,
				  WORLD_1D, &sig_send);
			// MPI-irecv -> 1
			MPI_Irecv(&bbuffs[(loop_step + 1) % 2].data[0],
				  bodies_per_rank * 4, MPI_DOUBLE, prev, 0,
				  WORLD_1D, &sig_recv);
			// compute
			lfi(bbuffs[loop_step % 2]);
			lfi.update();
		}

#ifdef VERBOSE
		objects.print_to_file(
		    "./out/bodies_r" + std::to_string(coord[0]), time_step);
#endif
		MPI_Wait(&sig_recv, NULL);
		MPI_Wait(&sig_send, NULL);
		MPI_Irecv(&bbuffs[0].data[0], bodies_per_rank * 4, MPI_DOUBLE,
			  prev, 0, WORLD_1D, &sig_recv);
		MPI_Isend(&objects.data[0], bodies_per_rank * 4, MPI_DOUBLE,
			  next, 0, WORLD_1D, &sig_send);
	}
	auto t_end = std::chrono::high_resolution_clock::now();
	auto t_diff = std::chrono::duration_cast<std::chrono::microseconds>(
			  t_end - t_start)
			  .count();

	if (coord[0] == 0) {
		std::cout << size << ", " << processes << ", " << t_diff << ", "
			  << iterations << "\n";
	}

	MPI::Finalize();
	return 0;
}

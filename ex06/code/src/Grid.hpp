#ifndef GRID_H
#define GRID_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct Grid {
	/*!
	 * @class simple helper struct for easy access via 2D coordinates while
	 * memory stays flat (no vector of vectors)
	 * _data is public so that it can be directly accessed (for MPI
	 * purposes)
	 */
	std::vector<double> _data;
	size_t dimx, dimy;

	// frobid default construction
	Grid() = delete;
	Grid(size_t dim);               // dim X dim grid
	Grid(size_t dimx, size_t dimy); // dimx X dimy grid
	virtual ~Grid();

	// use operator() since [] is painfull to implement for multi
	// dimensional access. NOTE: access pattern contrary to mathematical
	// notation!! idx is column index and idy row index to comply with
	// notation in exercise
	double &operator()(size_t idx, size_t idy);      // write
	double operator()(size_t idx, size_t idy) const; // read
	void to_csv(std::string filename, size_t iteration);
};

inline Grid::Grid(size_t dim) : dimx(dim), dimy(dim)
{
	_data.reserve(dimx * dimy);
	std::fill_n(std::back_inserter(_data), dimx * dimy, 0);
	size_t low = dimx / 4;
	size_t high = 3 * dimx / 4;
	for (size_t i = low; i < high; ++i) {
		_data[i] = 127.;
	}
}

inline Grid::Grid(size_t dimx, size_t dimy) : dimx(dimx), dimy(dimy)
{
	_data.reserve(dimx * dimy);
	std::fill_n(std::back_inserter(_data), dimx * dimy, 0);
	size_t low = dimx / 4;
	size_t high = 3 * dimx / 4;
	for (size_t i = low; i < high; ++i) {
		_data[i] = 127.;
	}
}

inline Grid::~Grid() {}

inline double &Grid::operator()(size_t idx, size_t idy)
{
	return _data[idy * dimx + idx];
}

inline double Grid::operator()(size_t idx, size_t idy) const
{
	return _data[idy * dimx + idx];
}

inline void Grid::to_csv(std::string filename, size_t iteration)
{
	std::ofstream fs;
	fs.open(filename, std::ios_base::app);
	bool first = true;
	fs << iteration << "; " << dimx << "; " << dimy << "; [";
	for (size_t j = 0; j < dimy; ++j) {
		for (size_t i = 0; i < dimx; ++i) {
			if (!first) {
				fs << ", ";
			}
			fs << _data[j * dimx + i];
			first = false;
		}
	}
	fs << "]\n";
	fs.flush();
}

#endif /* GRID_H */

#include <stdlib.h>
#include <math.h>
#include <string>
#include <valarray>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>

#include "export.hpp"
#include "timer.hpp"

void initializeRing(std::vector<double> & timed_grid, size_t size, uint32_t diameter, double heat = 127) {
    double cx, cy;
    cx = cy = size/2;

    for (uint32_t row =0; row < size; row++){
        for(uint32_t col = 0; col < size; col++){
            double d = sqrt(pow(abs(cx - col),2) + pow(abs(cy - row),2));
            if ( abs(d - diameter/2) < .5 ) timed_grid[col+row*size] = heat;
        }
    }
}

void iterate(std::vector<double> & t0_grid, std::vector<double> & t1_grid, size_t size) {
    for (uint32_t row = 1; row < size-1; row++){
        for(uint32_t col = 1; col < size-1; col++){

            t1_grid[col+row*size] = t0_grid[col+row*size]
                + 6./25 * (
                    -4 * t0_grid[col+row*size]
                    + t0_grid[col+1+row*size]
                    + t0_grid[col-1+row*size]
                    + t0_grid[col+(row+1)*size]
                    + t0_grid[col+(row-1)*size]
                );
        }
    }
}

int main(int argc, const char** argv) {
    std::vector<std::map<std::string,uint32_t>> givens {
        {{"size", 100}, {"diameter", 35}},
        {{"size", 500}, {"diameter", 175}},
        {{"size", 1000}, {"diameter", 350}},
        {{"size", 5000}, {"diameter", 1750}},
        {{"size", 10000}, {"diameter", 3500}},
    };

    constexpr size_t to_keep = 5;

    auto given = givens[atoi(argv[1])];
    std::cout << "running with size of "
        << given["size"] << "x" << given["size"] << "\n";
    std::array<std::vector<double>, to_keep> grid {}; // only a given no of timesteps will be kept

    for (auto & g: grid)
        g.resize(given["size"]*given["size"], 0); // col+row*given["size"]

    initializeRing(grid[0], given["size"], given["diameter"]);

    timer::Precision t1 {};
    for (uint32_t iter = 0; iter < 1000; iter++) {
        iterate(grid[iter%to_keep], grid[(iter+1)%to_keep], given["size"]);
        if (iter < 15){
            std::stringstream ss;
            ss << std::setfill('0') << "dat/" << given["size"] << "_t" << std::setw(4) << iter << ".pgm";
            Export::PGM pgm = Export::PGM(ss.str(), grid[iter%to_keep], given["size"], 127);
            pgm.flush();
        }

    }
    t1.stop();

    std::cout << "done with "
        << given["size"] << "x" << given["size"]
        << " (in " << t1.duration() << "s)\n";

    Export::Export ex = Export::Export("dat/dat.dat");
    ex.from_df({{given["size"]}, {t1.duration()}});
    ex.flush();
}

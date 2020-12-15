#include <chrono>
#include <iostream>
#include <random>
#include <string>
#include <array>
#include <vector>

#include "../../libs/argparse/argparse.hpp"
#include "../../libs/export/export.hpp"

int main(int argc, const char **argv) {
    // get args
    argparse::ArgumentParser ap;
    ap.addArgument("-s", "--size", 1, false);
    ap.parse(argc, argv);

    size_t size = ap.retrieve<size_t>("size");

    // setup buffers
    std::vector<std::vector<double>> buffers(2, std::vector<double>(size*size));

    Export::to_csv<double>("out/seq.csv", buffers, {"0", "1"});

    return 0;
}

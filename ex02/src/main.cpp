#include <string>
#include "../../libs/argparse/argparse.hpp"

int main(int argc,  const char** argv) {

    ArgumentParser ap;
    ap.parse(argc, argv);
    return 0;
}

#ifndef _EXPORT_HPP_
#define _EXPORT_HPP_

//#include <string>
#include <iostream>
//#include <iterator>
#include <fstream>
//#include <sstream>
#include <vector>
//#include <algorithm>
#include <array>

namespace Export{
    template <typename numerical_t, unsigned int col_n>
    void to_csv(std::string fp, std::vector<std::array<numerical_t, col_n>> &pixelarray, std::array<std::string, col_n> columnnames = {}) {
        std::ofstream fs;
        fs.open(fp);
        //header
        fs << "# ";
        const char *sep = "";
        for (auto column: columnnames) {
            fs << sep << column;
            sep = ", ";
        }
        //pixels
        for (auto row: pixelarray){
            fs << "\n";
            const char *sep = "";
            for (auto column: row) {
                fs << sep << column;
                sep = ", ";
            }
        }
        fs.close();
    }
}

#endif /* ifndef _EXPORT_HPP_ */

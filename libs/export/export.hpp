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
    template <typename numerical_t>
    void to_csv(
            std::string fp,
            std::vector<std::vector<numerical_t>> &columndata,
            std::vector<std::string> columnnames = {}) {

        std::cout << "exporting to " << fp << "\n";
        std::ofstream fs;
        fs.open(fp);
        //header
        fs << "# ";
        const char *sep = "";
        for (auto column: columnnames) {
            fs << sep << column;
            sep = ", ";
        }
        fs << "\n";

        size_t rows = 0;
        for (auto col: columndata){
            rows = (col.size() > rows) ? col.size() : rows;
        }

        for (size_t r = 0; r < rows; r++){
            const char *sep = "";
            for (auto col: columndata){
                //fs << "\n";
                fs << sep << col[r];
                sep = ", ";
            }
            fs << "\n";
        }
        fs.close();
    }
}

#endif /* ifndef _EXPORT_HPP_ */

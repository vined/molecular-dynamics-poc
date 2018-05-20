#include <fstream>
#include <iostream>

#include "Parameters.h"
#include "utils/FileUtils.h"

Parameters parseParams(std::vector<std::string> argv) {
    return Parameters(
            std::stof(argv[0]), // dt
            std::stof(argv[1]), // max time
            std::stof(argv[2]), // temperature
            std::stof(argv[3]), // sigma
            std::stof(argv[4]), // mass
            std::stof(argv[5]), // epsilon
            std::stol(argv[6]), // data export interval
            std::stof(argv[7]), // density
            std::stol(argv[8])  // lattice count
    );
}

bool isCommentOrEmpty(std::string line) {
    return line.size() == 0 || line[0] == '#';
}

std::string getValue(std::string line) {

    std::size_t pos = line.find_first_of(" #;,\t\n");

    if (pos != std::string::npos) {
        return line.substr(0, pos);
    }
    return line;
}

Parameters readParameters(std::string fileName) {

    std::ifstream paramsFile = openFile(fileName);
    std::string line;
    std::vector<std::string> params;

    while (getline(paramsFile, line)) {

        if (!isCommentOrEmpty(line)) {
            params.push_back(
                    getValue(line)
            );
        }
    }

    return parseParams(params);
}

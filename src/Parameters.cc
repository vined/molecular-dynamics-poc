#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "Parameters.h"

Parameters parseParams(std::vector<std::string> argv) {
    return Parameters(
            std::stof(argv[0]), // dt
            std::stof(argv[1]), // max time
            std::stof(argv[2]), // box size
            std::stol(argv[3]), // atoms count

            std::stof(argv[4]), // sigma
            std::stof(argv[5]), // mass
            std::stof(argv[6]), // epsilon
            std::stof(argv[7])  // energy cut-off
    );
}

std::ifstream openFile(std::string fileName) {
    std::ifstream file;
    file.open(fileName);

    if (file.is_open()) {
        return file;
    }

    std::cout << "Failed to open file: " << fileName << std::endl;
    throw 200;
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
    std::vector <std::string> params;

    while (getline(paramsFile, line)) {

        if (!isCommentOrEmpty(line)) {
            params.push_back(
                    getValue(line)
            );
        }
    }

    return parseParams(params);
}


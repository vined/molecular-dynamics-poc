#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "Parameters.h"

Parameters parseParams(std::vector<std::string> argv) {
    return Parameters(
            std::stof(argv[0]), // dt
            std::stof(argv[1]), // max time
            std::stof(argv[2]), // sigma
            std::stof(argv[3]), // mass
            std::stof(argv[4]), // epsilon
            std::stof(argv[5]), // energy cut-off
            std::stol(argv[6]), // atoms positions export interval
            std::stol(argv[7])  // atoms count
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

std::vector<std::string> splitBySpace(std::string atomData) {

    std::string delimiter = " \t";
    std::vector<std::string> data;

    std::size_t pos = 0;
    while ((pos = atomData.find_first_of(delimiter)) != std::string::npos) {
        data.push_back(atomData.substr(0, pos));
        atomData.erase(0, pos + 1);
    }
    data.push_back(atomData);
    return data;
}

Atom parseAtom(std::string atomData, AtomType atomType) {

    std::vector<std::string> params = splitBySpace(atomData);

    return Atom(atomType, {
            std::stof(params[1]) * pow(10, -9),
            std::stof(params[2]) * pow(10, -9),
            std::stof(params[3]) * pow(10, -9),
    });
}

std::vector<Atom> parseAtoms(std::vector<std::string> lines, AtomType atomType) {
    std::vector<Atom> atoms;

    for (std::string line : lines) {
        atoms.push_back(parseAtom(line, atomType));
    }

    return atoms;
}

std::vector<Atom> readAtomsData(Parameters params, std::string fileName) {

    std::ifstream paramsFile = openFile(fileName);

    std::string line;
    std::vector<std::string> lines;

    while (getline(paramsFile, line)) {
        lines.push_back(line);
    }

    AtomType atomType = {
            "Ar",
            params.mass,
            params.sigma,
            params.epsilon
    };

    return parseAtoms(lines, atomType);
}

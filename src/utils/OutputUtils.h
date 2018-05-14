#ifndef MOLECULAR_DYNAMICS_POC_OUTPUT_UTILS_H
#define MOLECULAR_DYNAMICS_POC_OUTPUT_UTILS_H

#include <vector>
#include <string>

#include "../Atoms.h"


std::string vectorToString(Vector v);

void exportVector(std::string name, std::vector<double> vect, int precision);

void exportMultiVector(
        std::string name,
        std::vector <std::string> colNames,
        std::vector <std::vector<double>> vectors,
        int precision
);

void exportAtomsPositions(std::string name, std::vector<Atom> atoms);

#endif //MOLECULAR_DYNAMICS_POC_OUTPUT_UTILS_H

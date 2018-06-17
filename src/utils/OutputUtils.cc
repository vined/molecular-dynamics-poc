#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>

#include "OutputUtils.h"
#include "../Molecules.h"


#define OUTPUT_DIR "./out/"
#define PRECISION 15

const double nmToAngst = 1.0;

std::string vectorToString(Vector v) {
    return std::to_string(v.x) + " " + std::to_string(v.y) + " " + std::to_string(v.z);
}

std::string quatToString(Quaternion q) {
    return std::to_string(q.a) + " " + std::to_string(q.b) + " " + std::to_string(q.c) + " " + std::to_string(q.d);
}

void exportVector(std::string name, std::vector<double> vect, int precision) {

    std::string fileName = OUTPUT_DIR + name + ".dat";

    std::cout << "Exporting " << fileName << std::endl;

    std::ofstream dat;
    dat.open(fileName);
    dat.precision(precision);

    dat << "# " << name << std::endl;

    for (double val : vect) {
        dat << val << std::endl;
    }

    dat.close();
    std::cout << "Exporting " << fileName << " done." << std::endl;
}

void exportMultiVector(
        std::string name,
        std::vector <std::string> colNames,
        std::vector <std::vector<double>> vectors,
        int precision
) {
    std::string fileName = OUTPUT_DIR + name + ".dat";

    std::cout << "Exporting vectors to " << fileName << std::endl;

    std::ofstream dat;
    dat.open(fileName);
    dat.precision(precision);

    dat << "# ";
    for (std::string val : colNames) {
        dat << val << '\t';
    }
    dat << std::endl;

    long n = vectors[0].size();
    for (std::vector<double> v: vectors) {
        if (v.size() < n) {
            n = v.size();
        }
    }

    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < vectors.size(); j++) {
            dat << vectors[j][i] << '\t';
        }
        dat << std::endl;
    }

    dat.close();
    std::cout << "Exporting vectors to " << fileName << " done." << std::endl;
}

void exportAtomsPositions(std::string name, std::vector<Atom> atoms) {
    std::string fileName = OUTPUT_DIR + name + ".xyz";

    std::ofstream dat;
    dat.open(fileName, std::ios_base::app);
    dat.precision(PRECISION);

    dat << atoms.size() << std::endl;

    dat << "# " << std::endl;

    for (Atom atom : atoms) {
        Vector p = atom.position;
        dat << atom.type.symbol << " " << p.x * nmToAngst << " " << p.y * nmToAngst << " " << p.z * nmToAngst << std::endl;
    }

    dat.close();
}

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

#include "Atoms.h"
#include "utils/FileUtils.h"

#define AR_MASS 6.69
#define AR_SIGMA 0.34
#define AR_EPSILON 1.65

AtomType argonAtomType = {
        "Ar",
        AR_MASS,
        AR_SIGMA,
        AR_EPSILON
};

std::vector<Atom> initializeAtoms(Parameters params, double *width) {

    std::vector<Atom> atoms;

    double half_step = params.density / 2.0;
    double half_half_step = half_step / 2.0;
    double x = half_half_step;

    for (int i = 0; i < params.lattice_count; i++) {

        double y = 0;
        for (int j = 0; j < params.lattice_count; j++) {

            double z = 0;
            for (int k = 0; k < params.lattice_count; k++) {
                atoms.push_back(Atom(argonAtomType, {x + half_step, y + half_step, z + half_step}));

//                atoms.push_back(Atom(atomType, {x, y, z}));
//                atoms.push_back(Atom(atomType, {x + half_half_step, y + half_half_step, z + half_half_step}));
//                atoms.push_back(Atom(atomType, {x, y + half_step, z + half_step}));
//                atoms.push_back(Atom(atomType, {x + half_step, y + half_step, z}));
//                atoms.push_back(Atom(atomType, {x + half_step, y, z + half_step}));
                z += params.density;
            }
            y += params.density;
        }
        x += params.density;
    }

    *width = x - half_half_step;
    return atoms;
}

std::vector<Atom> initializeVelocities(std::vector<Atom> atoms, double temperature) {

    Vector totalVelocity = getZeroVector();
    double velocityScale = sqrt(3.0 * (1.0 - (1.0 / atoms.size()) * temperature));

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::default_random_engine re;

    for (int i = 0; i < atoms.size(); i++) {
        atoms[i].velocity = scale({unif(re), unif(re), unif(re)}, velocityScale);
        totalVelocity = sum(totalVelocity, atoms[i].velocity);
    }

    Vector velocityReducer = scale(totalVelocity, -1.0 / (double) atoms.size());

    for (int i = 0; i < atoms.size(); i++) {
        atoms[i].velocity = sum(atoms[i].velocity, velocityReducer);
    }

    return atoms;
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
            std::stof(params[1]),
            std::stof(params[2]),
            std::stof(params[3]),
    });
}

std::vector<Atom> parseAtoms(std::vector<std::string> lines, AtomType atomType, double *box_size) {
    std::vector<Atom> atoms;

    for (std::string line : lines) {
        Atom a = parseAtom(line, atomType);
        atoms.push_back(a);

        double max = getMinMaxCoordinate(a.position).second;
        if (max > *box_size) {
            *box_size = max;
        }
    }

    return atoms;
}

std::vector<Atom> readAtomsData(Parameters params, std::string fileName, double *box_size) {

    std::ifstream paramsFile = openFile(fileName);

    std::string line;
    std::vector<std::string> lines;

    while (getline(paramsFile, line)) {
        lines.push_back(line);
    }

    return parseAtoms(lines, argonAtomType, box_size);
}

#ifndef MOLECULAR_DYNAMICS_POC_ATOMS_H
#define MOLECULAR_DYNAMICS_POC_ATOMS_H

#include <vector>
#include <string>

#include "Vectors.h"
#include "Parameters.h"


struct AtomType {
    std::string symbol;
    double mass;
    double sigma;
    double epsilon;
};

struct Atom {
    struct AtomType type;
    struct Vector position;
    struct Vector velocity;
    struct Vector acceleration;

    Atom(
            struct AtomType _type,
            struct Vector _position
    ) {
        type = _type;
        position = _position;
        velocity = getZeroVector();
        acceleration = getZeroVector();
    }
};

std::vector<Atom> initializeAtoms(Parameters params, double *width);

std::vector<Atom> initializeVelocities(std::vector<Atom> atoms, double velocityScale);

std::vector<Atom> readAtomsData(Parameters params, std::string fileName, double *box_size);

#endif //MOLECULAR_DYNAMICS_POC_ATOMS_H

#ifndef MOLECULAR_DYNAMICS_POC_ATOMS_H
#define MOLECULAR_DYNAMICS_POC_ATOMS_H

#include "Vectors.h"
#include "Parameters.h"


struct AtomType {
    int type = 0;
    double mass;
    double sigma;
    double epsilon;

    AtomType(
            double _mass,
            double _sigma,
            double _epsilon
    ) {
        mass = _mass;
        sigma = _sigma;
        epsilon = _epsilon;
    }
};

struct Atom {
    struct AtomType type;
    struct Vector position;
    struct Vector velocity;
    struct Vector potential;
    struct Vector acceleration;

    Atom(
            struct AtomType _type,
            struct Vector _position
    ) {
        type = _type;
        position = _position;
        velocity = getZeroVector();
        potential = getZeroVector();
        acceleration = getZeroVector();
    }
};

std::vector<Atom> initializeAtoms(Parameters params);

#endif //MOLECULAR_DYNAMICS_POC_ATOMS_H

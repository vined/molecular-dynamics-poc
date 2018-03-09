#ifndef MOLECULAR_DYNAMICS_POC_ATOMS_H
#define MOLECULAR_DYNAMICS_POC_ATOMS_H

#include "vectors.h"


struct AtomType {
    int type;
    double mass;
};

struct Atom {
    struct AtomType type;
    struct Vector position;
    struct Vector velocity;
    struct Vector potential;
    struct Vector acceleration;
};


#endif //MOLECULAR_DYNAMICS_POC_ATOMS_H



#ifndef MOLECULAR_DYNAMICS_POC_MOLECULE_H
#define MOLECULAR_DYNAMICS_POC_MOLECULE_H

#include <vector>

#include "Atoms.h"
#include "Sites.h"
#include "Vectors.h"
#include "Quaternions.h"
#include "Matrix.h"


struct Molecule {
    std::vector<Site> sites;
    Vector acceleration = getZeroVector();
    Vector position = getZeroVector();
    Vector torque = getZeroVector();
    Vector velocity = getZeroVector();
    Quaternion quaternion = getZero();
    Quaternion qVelocity = getZero();
    Quaternion qAcceleration = getZero();
};

Vector computeAngularVelocities(Molecule m);


#endif //MOLECULAR_DYNAMICS_POC_MOLECULE_H

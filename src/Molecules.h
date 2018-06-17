

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
    Vector acceleration1 = getZeroVector();
    Vector acceleration2 = getZeroVector();
    Vector position = getZeroVector();
    Vector positionPrev = getZeroVector();
    Vector torque = getZeroVector();
    Vector velocity = getZeroVector();
    Vector velocityPrev = getZeroVector();
    Vector inertia = getZeroVector();
    Quaternion quaternion = getZero();
    Quaternion quaternionPrev = getZero();
    Quaternion qVelocity = getZero();
    Quaternion qVelocityPrev = getZero();
    Quaternion qAcceleration = getZero();
    Quaternion qAcceleration1 = getZero();
    Quaternion qAcceleration2 = getZero();
};

Vector computeAngularVelocities(Molecule m);

void computeAccelerationQuats(std::vector<Molecule> *molecules);

void calculateTorques(std::vector<Molecule> *molecules);

void updateSitesCoordinates(std::vector<Molecule> *molecules);

void adjustQuaternions(std::vector<Molecule> *molecules);

#endif //MOLECULAR_DYNAMICS_POC_MOLECULE_H

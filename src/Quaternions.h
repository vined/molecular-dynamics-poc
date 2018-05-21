#ifndef MOLECULAR_DYNAMICS_POC_QUATERNION_H
#define MOLECULAR_DYNAMICS_POC_QUATERNION_H

#include <vector>


struct Quaternion {
    double a, b, c, d;
};


Quaternion getZero();

Quaternion add(Quaternion q1, Quaternion q2);

Quaternion scale(Quaternion q, double scalar);

Quaternion multiply(Quaternion q1, Quaternion q2);

double length(Quaternion q);

double squareLength(Quaternion q);

std::vector<double> toVector(Quaternion q);

#endif //MOLECULAR_DYNAMICS_POC_QUATERNION_H

#ifndef MOLECULAR_DYNAMICS_POC_QUATERNION_H
#define MOLECULAR_DYNAMICS_POC_QUATERNION_H


struct Quaternion {
    double a, b, c, d;
};


Quaternion getZero();

Quaternion add(Quaternion q1, Quaternion q2);

Quaternion scale(Quaternion q, double scalar);

Quaternion multiply(Quaternion q1, Quaternion q2);

double length(Quaternion q);

#endif //MOLECULAR_DYNAMICS_POC_QUATERNION_H

#ifndef MOLECULAR_DYNAMICS_POC_VECTOR_UTILS_H
#define MOLECULAR_DYNAMICS_POC_VECTOR_UTILS_H


struct Vector {
    double x;
    double y;
    double z;
};

Vector getZeroVector();

Vector subtract(Vector v1, Vector v2);

Vector sum(Vector v1, Vector v2);

Vector sum(Vector v1, Vector v2, Vector v3);

Vector add(Vector v, double d);

Vector scale(Vector v, double scalar);

Vector multiply(Vector v1, Vector v2);

double squaredLength(Vector v);

double length(Vector v);

std::pair<double, double> getMinMaxCoordinate(Vector v);


#endif //MOLECULAR_DYNAMICS_POC_VECTOR_UTILS_H

#ifndef MOLECULAR_DYNAMICS_POC_VECTOR_UTILS_H
#define MOLECULAR_DYNAMICS_POC_VECTOR_UTILS_H


struct Vector {
    double x;
    double y;
    double z;
};

struct Vector getZeroVector();

struct Vector subtract(struct Vector v1, struct Vector v2);

struct Vector sum(struct Vector v1, struct Vector v2);

struct Vector sum(struct Vector v1, struct Vector v2, struct Vector v3);

struct Vector add(struct Vector v, double d);

struct Vector scale(struct Vector v, double scalar);

double squaredLength(struct Vector v);

double length(struct Vector v);


#endif //MOLECULAR_DYNAMICS_POC_VECTOR_UTILS_H

#include <cmath>

#include "Vectors.h"


// TODO add unit tests

struct Vector subtract(struct Vector v1, struct Vector v2) {
    return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

struct Vector sum(struct Vector v1, struct Vector v2) {
    return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

struct Vector sum(struct Vector v1, struct Vector v2, struct Vector v3) {
    return {v1.x + v2.x + v3.x, v1.y + v2.y + v3.y, v1.z + v2.z + v3.z};
}

struct Vector add(struct Vector v, double d) {
    return {v.z + d, v.y + d, v.z + d};
}

struct Vector scale(struct Vector v, double scalar) {
    return {v.x * scalar, v.y * scalar, v.z * scalar};
}

double squaredLength(struct Vector v) {
    return std::pow(v.x, 2.0) + std::pow(v.y, 2.0) + std::pow(v.z, 2.0);
}

double length(struct Vector v) {
    return sqrt(std::pow(v.x, 2.0) + std::pow(v.y, 2.0) + std::pow(v.z, 2.0));
}

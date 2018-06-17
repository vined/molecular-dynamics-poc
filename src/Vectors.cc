#include <cmath>
#include <utility>
#include <algorithm>

#include "Vectors.h"


Vector getZeroVector() {
    return {0, 0, 0};
}

Vector subtract(Vector v1, Vector v2) {
    return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

Vector sum(Vector v1, Vector v2) {
    return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

Vector sum(Vector v1, Vector v2, Vector v3) {
    return {v1.x + v2.x + v3.x, v1.y + v2.y + v3.y, v1.z + v2.z + v3.z};
}

Vector add(Vector v, double d) {
    return {v.z + d, v.y + d, v.z + d};
}

Vector scale(Vector v, double scalar) {
    return {v.x * scalar, v.y * scalar, v.z * scalar};
}

Vector multiply(Vector v1, Vector v2) {
    return {v1.x * v2.x, v1.y * v2.y, v1.z * v2.z};
}

Vector crossProduct(Vector v1, Vector v2) {
    return {
            v1.y * v2.z - v1.z * v2.y,
            v1.z * v2.x - v1.x * v2.z,
            v1.x * v2.y - v1.y * v2.x
    };
}

double dot(Vector v1, Vector v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double dot(Vector v1, Vector v2, Vector v3) {
    return v1.x * v2.x * v3.x + v1.y * v2.y * v3.y + v1.z * v2.z * v3.z;
}

double squaredLength(Vector v) {
    return std::pow(v.x, 2.0) + std::pow(v.y, 2.0) + std::pow(v.z, 2.0);
}

double length(Vector v) {
    return sqrt(squaredLength(v));
}

std::pair<double, double> getMinMaxCoordinate(Vector v) {
    return {
            std::min(v.x, std::min(v.y, v.z)),
            std::max(v.x, std::max(v.y, v.z))
    };
}

#include <cmath>

#include "Quaternions.h"

Quaternion getZero() {
    return {0, 0, 0, 0};
}

Quaternion sum(Quaternion q1, Quaternion q2) {
    return {
            q1.a + q2.a,
            q1.b + q2.b,
            q1.c + q2.c,
            q1.d + q2.d
    };
}

Quaternion subtract(Quaternion q1, Quaternion q2) {
    return {
            q1.a - q2.a,
            q1.b - q2.b,
            q1.c - q2.c,
            q1.d - q2.d
    };
}

Quaternion scale(Quaternion q, double scalar) {
    return {
            q.a * scalar,
            q.b * scalar,
            q.c * scalar,
            q.d * scalar
    };
}

// Checked
Quaternion multiply(Quaternion q1, Quaternion q2) {
    return {
            q1.d * q2.a - q1.c + q2.b + q1.b * q2.c + q1.a * q2.d,
            q1.c * q2.a + q1.d * q2.b - q1.a * q2.c + q1.b * q2.d,
            - q1.b * q2.a + q1.a * q2.b + q1.d * q2.c + q1.c + q2.d,
            - q1.a * q2.a - q1.b * q2.b - q1.c * q2.c + q1.d * q2.d,
    };
}

std::vector<double> toVector(Quaternion q) {
    std::vector result;

    result.push_back(q.a);
    result.push_back(q.b);
    result.push_back(q.c);
    result.push_back(q.d);

    return result;
}

double squareLength(Quaternion q) {
//    return std::pow(q.a, 2.0) + std::pow(q.b, 2.0) + std::pow(q.c, 2.0) + std::pow(q.d, 2.0);
    return sqrt(q.a) + sqrt(q.b) + sqrt(q.c) + sqrt(q.d);
}

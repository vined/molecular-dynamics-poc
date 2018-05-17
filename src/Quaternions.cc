#include <cmath>

#include "Quaternions.h"

Quaternion getZero() {
    return {0, 0, 0, 0};
}

Quaternion add(Quaternion q1, Quaternion q2) {
    return {
            q1.a + q2.a,
            q1.b + q2.b,
            q1.c + q2.c,
            q1.d + q2.d
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

Quaternion multiply(Quaternion q1, Quaternion q2) {
    return {
            q1.a * q2.a - q1.b * q2.b - q1.c * q2.c - q1.d * q2.d,
            q1.a * q2.b + q1.b * q2.a - q1.c * q2.d + q1.d * q2.c,
            q1.a * q2.c + q1.b * q2.d + q1.c * q2.a - q1.d * q2.b,
            q1.a * q2.d - q1.b * q2.c + q1.c * q2.b + q1.d * q2.a,
    };
}

double length(Quaternion q) {
    return sqrt(std::pow(q.a, 2.0) + std::pow(q.b, 2.0) + std::pow(q.c, 2.0) + std::pow(q.d, 2.0));
}
#include "Matrix.h"


Vector multiplyVector(Matrix mat, Vector v) {
    std::vector<double> m = mat.values;

    return {
        m[0] * v.x + m[3] * v.y + m[6] * v.z,
        m[1] * v.x + m[4] * v.y + m[7] * v.z,
        m[2] * v.x + m[5] * v.y + m[8] * v.z,
    };
}

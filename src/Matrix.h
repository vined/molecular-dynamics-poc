#ifndef MOLECULAR_DYNAMICS_POC_MATRIX_H
#define MOLECULAR_DYNAMICS_POC_MATRIX_H

#include <vector>

#include "Vectors.h"


struct Matrix {
    // Should be size of 9
    std::vector<double> values;
};

Vector multiplyVector(Matrix m, Vector v);


#endif //MOLECULAR_DYNAMICS_POC_MATRIX_H

#ifndef MOLECULAR_DYNAMICS_POC_PARAMETERS_H
#define MOLECULAR_DYNAMICS_POC_PARAMETERS_H

#include <vector>

#include "atoms.h"
#include "vectors.h"


struct Parameters {
    int dt; // delta t - time step size (~10e-15s)
    double max_time; // maximal simulation time (1-10ns

    double model_size_x;
    double model_size_y;
    double model_size_z;

    int atomsCnt; // atoms count in simulation
    std::vector<Atom> atoms;
};


#endif //MOLECULAR_DYNAMICS_POC_PARAMETERS_H

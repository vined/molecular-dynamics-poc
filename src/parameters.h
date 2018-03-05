

#ifndef MOLECULAR_DYNAMICS_POC_PARAMETERS_H
#define MOLECULAR_DYNAMICS_POC_PARAMETERS_H

struct Parameters {
    int dt; // delta t - time step size (10^-15s)
    double maxTime; // maximal simulation time
    int atomsCnt; // atoms count in simulation
};

#endif //MOLECULAR_DYNAMICS_POC_PARAMETERS_H

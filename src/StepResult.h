

#ifndef MOLECULAR_DYNAMICS_POC_STEPRESULT_H
#define MOLECULAR_DYNAMICS_POC_STEPRESULT_H

#include "Atoms.h"

struct StepResult {
    double potentialEnergy;
    // Todo probably remove
    double force;
    Molecule m1;
    Molecule m2;
};

#endif //MOLECULAR_DYNAMICS_POC_STEPRESULT_H

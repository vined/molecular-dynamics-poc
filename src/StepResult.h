

#ifndef MOLECULAR_DYNAMICS_POC_STEPRESULT_H
#define MOLECULAR_DYNAMICS_POC_STEPRESULT_H

#include "Atoms.h"

struct StepResult {
    double potentialEnergy;
    // Todo probably remove
    double force;
};

#endif //MOLECULAR_DYNAMICS_POC_STEPRESULT_H

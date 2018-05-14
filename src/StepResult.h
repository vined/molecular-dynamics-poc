

#ifndef MOLECULAR_DYNAMICS_POC_STEPRESULT_H
#define MOLECULAR_DYNAMICS_POC_STEPRESULT_H

#include "Atoms.h"

struct StepResult {
    double potentialEnergy;
    Atom a1;
    Atom a2;
};

#endif //MOLECULAR_DYNAMICS_POC_STEPRESULT_H

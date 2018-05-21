

#ifndef MOLECULAR_DYNAMICS_POC_SITES_H
#define MOLECULAR_DYNAMICS_POC_SITES_H

#include <string>

#include "Vectors.h"


struct Site{
    int type;
    Vector force;
    Vector position;
};


#endif //MOLECULAR_DYNAMICS_POC_SITES_H

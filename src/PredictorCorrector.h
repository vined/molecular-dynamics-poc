

#ifndef MOLECULAR_DYNAMICS_POC_PREDICTORCORRECTOR_H
#define MOLECULAR_DYNAMICS_POC_PREDICTORCORRECTOR_H

#include "Molecules.h"
#include "Quaternions.h"


void predict(std::vector<Molecule> *molecules, double dt);
void predictQuaternion(std::vector<Molecule> *molecules, double dt);

void correct(std::vector<Molecule> *molecules, double dt);
void correctQuaternion(std::vector<Molecule> *molecules, double dt);


#endif //MOLECULAR_DYNAMICS_POC_PREDICTORCORRECTOR_H

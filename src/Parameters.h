#ifndef MOLECULAR_DYNAMICS_POC_PARAMETERS_H
#define MOLECULAR_DYNAMICS_POC_PARAMETERS_H

#include <cmath>
#include <vector>
#include <string>

#include "Atoms.h"
#include "Vectors.h"


struct Parameters {

    double dt; // delta t - time step size in femtoseconds (e-15)
    double max_time; // maximum simulation time in nanoseconds (e-9)

    double box_size; //Model size in nanometers (e-9) (~31 A*)
    long atoms_cnt; // atoms count in simulation

    double sigma; // sigma in nano ergs (e-9)
    double mass; // atom mass in e-23 grams
    double epsilon; // epsilon in e-14
    double cutoff; // energy cut-off
    double unit_size; // temporary: atoms unit box size


    Parameters(
            double _dt, // femtoseconds (e-15)
            double _max_time, // nanoseconds (e-9)
            double _box_size, //nanometers (e-9)
            long _atoms_cnt,

            double _sigma, // nano ergs (e-9)
            double _mass, // e-23 grams
            double _epsilon, // e-14
            double _cutoff,
            double _unit_size
    ) {
        dt = _dt * pow(10, -15);
        max_time = _max_time * pow(10, -9);
        box_size = _box_size * pow(10, -9);
        atoms_cnt = _atoms_cnt;

        sigma = _sigma * pow(10, -9);
        mass = _mass * pow(10, -23);
        epsilon = _epsilon * pow(10, -14);
        cutoff = _cutoff;
        unit_size = _unit_size * pow(10, -9);
    }

};

Parameters readParameters(std::string fileName);

#endif //MOLECULAR_DYNAMICS_POC_PARAMETERS_H

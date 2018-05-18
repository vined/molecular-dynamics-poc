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

    double sigma; // sigma in nano ergs (e-9)
    double mass; // atom mass in e-23 grams
    double epsilon; // epsilon in e-14
    double cutoff; // energy cut-off
    long pos_export_interval; // atoms positions export interval
    long atoms_count; // atoms positions export interval


    Parameters(
            double _dt, // femtoseconds (e-15)
            double _max_time, // nanoseconds (e-9)
            double _sigma, // nano ergs (e-9)
            double _mass, // e-26 grams
            double _epsilon, // e-21
            double _cutoff, //nanometers (e-9)
            long _pos_export_interval,
            long _atoms_count
    ) {
        dt = _dt * pow(10, -15);
        max_time = _max_time * pow(10, -9);
        sigma = _sigma * pow(10, -9);
        mass = _mass * pow(10, -26);
        epsilon = _epsilon * pow(10, -21);
        cutoff = _cutoff * pow(10, -9);
        pos_export_interval = _pos_export_interval;
        atoms_count = _atoms_count;
    }

};

Parameters readParameters(std::string fileName);

std::vector<Atom> readAtomsData(Parameters params, std::string fileName);

#endif //MOLECULAR_DYNAMICS_POC_PARAMETERS_H

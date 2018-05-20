#ifndef MOLECULAR_DYNAMICS_POC_PARAMETERS_H
#define MOLECULAR_DYNAMICS_POC_PARAMETERS_H

#include <cmath>
#include <vector>
#include <string>

#include "Vectors.h"


struct Parameters {

    double dt; // delta t - time step size in femtoseconds
    double max_time; // maximum simulation time in nanoseconds
    double temperature;
    double sigma; // sigma in nano ergs (e-9)
    double mass; // atom mass in e-23 grams
    double epsilon; // epsilon in e-14
    long data_export_interval; // atoms positions export interval
    double density; // density
    long lattice_count;


    Parameters(
            double _dt,
            double _max_time,
            double _temperature,
            double _sigma, // nano ergs (e-9)
            double _mass, // e-26 grams
            double _epsilon, // e-21
            long _pos_export_interval,
            double _density,
            long _lattice_count
    ) {
        dt = _dt;
        max_time = _max_time;
        temperature =
        sigma = _sigma * pow(10, -9);
        mass = _mass * pow(10, -26);
        epsilon = _epsilon * pow(10, -21);
        data_export_interval = _pos_export_interval;
        density = _density;
        lattice_count = _lattice_count;
    }

};

Parameters readParameters(std::string fileName);

#endif //MOLECULAR_DYNAMICS_POC_PARAMETERS_H

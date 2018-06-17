#ifndef MOLECULAR_DYNAMICS_POC_PARAMETERS_H
#define MOLECULAR_DYNAMICS_POC_PARAMETERS_H

#include <cmath>
#include <vector>
#include <string>

#include "Vectors.h"


struct Parameters {

    double dt;
    double max_time;

    double temperature;
    double cut_off; // epsilon in e-14
    double density; // density
    long lattice_count;

    long data_export_interval;
    long adjust_temperature_interval;
    long adjust_equilibration_temp_interval;
    long equilibration_steps; // how long equilibration should be applied

    Parameters(
            double _dt,
            double _max_time,
            double _temperature,
            double _cut_off, // nano ergs (e-9)
            double _density,
            long _lattice_count,
            long _pos_export_interval,
            long _adjust_temperature_interval,
            long _adjust_equilibration_temp_interval,
            long _equilibration_steps
    ) {
        dt = _dt;
        max_time = _max_time;

        temperature = _temperature;
        cut_off = _cut_off;
        density = _density;
        lattice_count = _lattice_count;

        data_export_interval = _pos_export_interval;
        adjust_temperature_interval = _adjust_temperature_interval;
        adjust_equilibration_temp_interval = _adjust_equilibration_temp_interval;
        equilibration_steps = _equilibration_steps;
    }

};

Parameters readParameters(std::string fileName);

#endif //MOLECULAR_DYNAMICS_POC_PARAMETERS_H

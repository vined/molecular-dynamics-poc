

#ifndef MOLECULAR_DYNAMICS_POC_ENERGIES_H
#define MOLECULAR_DYNAMICS_POC_ENERGIES_H

struct Energies {
    double kinetic;
    double potential;
    double pressure;

    Energies(
            double _kinetic,
            double _potential,
            double _pressure
    ) {
        kinetic = _kinetic;
        potential = _potential;
        pressure = _pressure;
    }
};

#endif //MOLECULAR_DYNAMICS_POC_ENERGIES_H

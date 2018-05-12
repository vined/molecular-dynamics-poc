

#ifndef MOLECULAR_DYNAMICS_POC_ENERGIES_H
#define MOLECULAR_DYNAMICS_POC_ENERGIES_H

struct Energies {
    double kinetic;
    double potential;

    Energies(
            double _kinetic,
            double _potential
    ) {
        kinetic = _kinetic;
        potential = _potential;
    }
};

#endif //MOLECULAR_DYNAMICS_POC_ENERGIES_H

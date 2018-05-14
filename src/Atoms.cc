#include <cmath>

#include "Atoms.h"


std::vector<Atom> initializeAtoms(Parameters params) {

    AtomType atomType = {
            "Ar",
            params.mass,
            params.sigma,
            params.epsilon
    };

    std::vector<Atom> atoms;

    double x = 0;
    double half = params.unit_size / 2;

    while (x < params.box_size - half) {

        double y = 0;
        while (y < params.box_size - half) {

            double z = 0;
            while (z < params.box_size - half) {
                z += params.unit_size;
                atoms.push_back(Atom(atomType, {x, y, z}));
                atoms.push_back(Atom(atomType, {x + half, y + half, z + half}));
//                atoms.push_back(Atom(atomType, {x, y, z}));
//                atoms.push_back(Atom(atomType, {x, y + half, z + half}));
//                atoms.push_back(Atom(atomType, {x + half, y + half, z}));
//                atoms.push_back(Atom(atomType, {x + half, y, z + half}));
            }
            y += params.unit_size;
        }
        x += params.unit_size;
    }

    return atoms;
}

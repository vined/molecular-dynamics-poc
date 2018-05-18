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
    double border_length = 1.218 * params.atoms_count;
    double step_size = border_length / sqrt(params.atoms_count);
    double half_step = step_size / 2;

    while (x < border_length - half_step) {

        double y = 0;
        while (y < border_length - half_step) {

            double z = 0;
            while (z < border_length - half_step) {
                z += step_size;
                atoms.push_back(Atom(atomType, {x, y, z}));
                atoms.push_back(Atom(atomType, {x + half_step, y + half_step, z + half_step}));
//                atoms.push_back(Atom(atomType, {x, y, z}));
//                atoms.push_back(Atom(atomType, {x, y + half, z + half}));
//                atoms.push_back(Atom(atomType, {x + half, y + half, z}));
//                atoms.push_back(Atom(atomType, {x + half, y, z + half}));
            }
            y += step_size;
        }
        x += step_size;
    }

    return atoms;
}

#include <cmath>
#include <iostream>

#include "Atoms.h"
#include "Energies.h"
#include "Parameters.h"
#include "StepResult.h"
#include "utils/OutputUtils.h"

// Boltzmann constant erg*K
const double K_B = 1.38062e-16;


double checkBoundary(double position, double box_size) {
    if (position < 0) {
        return position - box_size;
    } else if (position > box_size) {
        return position - box_size;
    }

    return position;
}

Vector checkBoundaries(Vector position, double box_size) {
    return {
            checkBoundary(position.x, box_size),
            checkBoundary(position.y, box_size),
            checkBoundary(position.z, box_size)
    };
}

double withBoundary(double pos, double box_size, double cutoff) {
    return pos > cutoff ? pos : pos + box_size;
}

double getAxialDistance(double pos1, double pos2, double box_size, double cutoff) {
    return withBoundary(pos1, box_size, cutoff) - withBoundary(pos2, box_size, cutoff);
}

Vector getDistance(Vector pos1, Vector pos2, double box_size, double cutoff) {

    return {
            getAxialDistance(pos1.x, pos2.x, box_size, cutoff),
            getAxialDistance(pos1.y, pos2.y, box_size, cutoff),
            getAxialDistance(pos1.z, pos2.z, box_size, cutoff)
    };
}

//Lenard-Jones Potential for van der Waals system
StepResult getPotentialAndUpdateForEach(Atom a1, Atom a2, double box_size, double cutoff) {

    double sigma_squared = pow(a1.type.sigma, 2.0),
            cutoff_squared = pow(cutoff, 2.0),
            repulsion_erg = 48.0 * a1.type.epsilon,
            system_repulsion_erg = 4.0 * a1.type.epsilon;

    double result_potential = 0;

    Vector diff = getDistance(a1.position, a2.position, box_size, cutoff);
    double length_squared = std::pow(diff.x, 2.0) + std::pow(diff.y, 2.0) + std::pow(diff.z, 2.0);

    std::cout.precision(15);

    if (length_squared < cutoff_squared) {

        // Temporary variables to speed up calculations
        double a_squared = sigma_squared / length_squared, // (sigma / r_ij) ^ 2
                a_sixth = std::pow(a_squared, 3.0); // (sigma / r_ij) ^ 2

        double potential = repulsion_erg * a_sixth * (a_sixth - 0.5) / length_squared;
        Vector pairPotential = scale(diff, potential);

        // Update atoms potential energies
        a1.potential = sum(a1.potential, pairPotential);
        a2.potential = subtract(a2.potential, pairPotential);

        // return system potential
        result_potential = system_repulsion_erg * a_sixth * (a_sixth - 1.0);
    }
    return {result_potential, a1, a2};
}

double getTotalPotentialAndUpdateForEach(std::vector<Atom> *atoms, double box_size, double cutoff) {

    double potential_energy = 0.0;

    for (int i = 0; i < (*atoms).size(); i++) {
        for (int j = i + 1; j < (*atoms).size(); j++) {
            StepResult sr = getPotentialAndUpdateForEach((*atoms)[i], (*atoms)[j], box_size, cutoff);
            potential_energy += sr.potentialEnergy;
            (*atoms)[i] = sr.a1;
            (*atoms)[j] = sr.a2;
        }
    }

    return potential_energy;
}

// Velocity-Verlet method
Energies velocityVerlet(std::vector<Atom> *atoms, double dt, double box_size, double cutoff) {

    double mass = (*atoms)[0].type.mass,
            n = (*atoms).size(),
            half_dt = dt / 2.0,
            mass_inv = 1 / mass;

    // todo:
    // - check if atom is part of neigbouring chain
    // - eval forces of bonded atoms

    // calculate velocity and position for half step
    for (int i = 0; i < n; i++) {
        Atom a = (*atoms)[i];
        a.velocity = sum(a.velocity, scale(a.acceleration, half_dt));
        a.position = checkBoundaries(sum(a.position, scale(a.velocity, dt)), box_size);
        (*atoms)[i] = a;
    }

    double potential_energy = getTotalPotentialAndUpdateForEach(atoms, box_size, cutoff),
            kinetic_energy = 0;

    for (int i = 0; i < n; i++) {
        Atom a = (*atoms)[i];
        a.acceleration = scale(a.potential, mass_inv);
        a.velocity = sum(a.velocity, scale(a.acceleration, half_dt));
        (*atoms)[i] = a;
        kinetic_energy += squaredLength(a.velocity);
    }

    kinetic_energy = mass * kinetic_energy / 2.0;

    return Energies(kinetic_energy, potential_energy);
}

double getTemperature(double kineticEnergy, long count) {
    return (2.0 * kineticEnergy) / (3 * K_B * (double) count);
}

void runSimulation(Parameters params, std::vector<Atom> atoms) {

    double time = 0;
    long i = 0;

    std::vector<std::vector<double>> energies;
    std::vector<double> kinetic;
    std::vector<double> potential;
    energies.push_back(kinetic);
    energies.push_back(potential);

    std::vector<double> temperatures;
    std::vector<Vector> positions;

    while (time < params.max_time) {

        Energies e = velocityVerlet(&atoms, params.dt, params.box_size, params.cutoff);

        energies[0].push_back(e.kinetic);
        energies[1].push_back(e.potential);
        // TODO: use params.atoms_count
        temperatures.push_back(getTemperature(e.kinetic, atoms.size()));

        if (i % params.pos_export_interval  == 0) {
            exportAtomsPositions("atoms", atoms);
            std::cout << ".";
        }

        time += params.dt;
        i++;
    }

    std::cout << std::endl;

    std::cout << "Exporting results" << std::endl;
    std::vector<std::string> energiesCols;
    energiesCols.push_back("kinetic");
    energiesCols.push_back("potential");

    exportMultiVector("energies", energiesCols, energies, 16);
    exportVector("temperature", temperatures, 16);
}

int main(int argc, char *argv[]) {

    std::cout.precision(15);
    std::cout << "Reading parameters..." << std::endl;
    Parameters params = readParameters(argv[1]);

    std::vector<Atom> atoms;
    if (argc > 2) {
        std::cout << "Reading atoms..." << std::endl;
        atoms = readAtomsData(params, argv[2]);
    } else {
        std::cout << "Initializing atoms..." << std::endl;
        atoms = initializeAtoms(params);
    }

    std::cout << "Initialized " << atoms.size() << " atoms" << std::endl;

    std::cout << "Running simulation" << std::endl;
    runSimulation(params, atoms);
    
    return 0;
}



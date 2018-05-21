#include <cmath>
#include <iostream>
#include <utility>

#include "Atoms.h"
#include "Energies.h"
#include "Parameters.h"
#include "StepResult.h"
#include "LeapFrogResult.h"
#include "utils/OutputUtils.h"


double applyBoundary(double position, double box_size) {
    if (position < 0) {
        return position + box_size;
    } else if (position >= box_size) {
        return position - box_size;
    }

    return position;
}

Vector applyBoundaries(Vector position, double box_size) {
    return {
            applyBoundary(position.x, box_size),
            applyBoundary(position.y, box_size),
            applyBoundary(position.z, box_size)
    };
}

double withBoundary(double pos, double box_size, double cutoff) {
    return pos >= cutoff ? pos : pos + box_size;
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

    double cutoff_squared = pow(cutoff, 2.0);
    double result_potential = 0;
    double force = 0;

    Vector diff = getDistance(a1.position, a2.position, box_size, cutoff);
    double length_squared = std::pow(diff.x, 2.0) + std::pow(diff.y, 2.0) + std::pow(diff.z, 2.0);

    std::cout.precision(15);

    if (length_squared < cutoff_squared) {

        // Temporary variable to speed up calculations
        double len_inv_cubed = std::pow(1.0 / length_squared, 3.0);

        force = (48.0 * len_inv_cubed * (len_inv_cubed - 0.5)) / length_squared;
        Vector pairPotential = scale(diff, force);

        // Update atoms potential energies
        a1.acceleration = sum(a1.acceleration, pairPotential);
        a2.acceleration = subtract(a2.acceleration, pairPotential);

        // return system potential
        result_potential = 4.0 * len_inv_cubed * (len_inv_cubed - 1.0) + 1.0;
    }
    return {result_potential, force * length_squared, a1, a2};
}

LeapFrogResult getTotalPotentialAndUpdateForEach(std::vector<Atom> *atoms, double box_size, double cutoff) {

    double potential_energy = 0.0;
    double forces = 0.0;

    for (int i = 0; i < (*atoms).size(); i++) {
        for (int j = i + 1; j < (*atoms).size(); j++) {
            StepResult sr = getPotentialAndUpdateForEach((*atoms)[i], (*atoms)[j], box_size, cutoff);
            potential_energy += sr.potentialEnergy;
            forces += sr.force;
            (*atoms)[i] = sr.a1;
            (*atoms)[j] = sr.a2;
        }
    }

    return {potential_energy, forces};
}

Energies leapFrog(std::vector<Atom> *atoms, double dt, double box_size, double cutoff, double density) {

    double n = (*atoms).size();
    double half_dt = dt / 2.0;

    // todo:
    // - check if atom is part of neigbouring chain
    // - eval forces of bonded atoms

    for (int i = 0; i < n; i++) {
        Atom a = (*atoms)[i];
        a.velocity = sum(a.velocity, scale(a.acceleration, half_dt));
        a.position = applyBoundaries(sum(a.position, scale(a.velocity, dt)), box_size);
        std::pair<double, double> min_max= getMinMaxCoordinate(a.position);

        if (min_max.first < 0 || min_max.second > box_size) {
            Vector p = a.position;
            std::cout << "\nAtom " << i << " out of bounds (" << p.x << ',' << p.y << ',' << p.z << "), box size: " << box_size << std::endl;
            throw 100;
        }

        a.acceleration = getZeroVector();
        (*atoms)[i] = a;
    }

    // calculate velocity and position for half step
    LeapFrogResult lfResult = getTotalPotentialAndUpdateForEach(atoms, box_size, cutoff);

    double kinetic_energy = 0;
    for (int i = 0; i < n; i++) {
        Atom a = (*atoms)[i];
        a.velocity = sum(a.velocity, scale(a.acceleration, half_dt));
        (*atoms)[i] = a;
        kinetic_energy += squaredLength(a.velocity);
    }

    double pressure = density * (kinetic_energy + lfResult.forces) / ((*atoms).size() * 3.0);

    return Energies(
            kinetic_energy / 2.0,
            lfResult.totalPotential,
            pressure
    );
}

void runSimulation(Parameters params, std::vector<Atom> atoms, double box_size, double cutoff) {

    double time = 0.0;
    long i = 0;

    std::vector<std::vector<double>> energies;
    std::vector<double> kinetic;
    std::vector<double> potential;
    std::vector<double> pressure;
    energies.push_back(kinetic);
    energies.push_back(potential);
    energies.push_back(pressure);

    std::vector<double> temperatures;
    std::vector<Vector> positions;

    while (time < params.max_time) {

        Energies e = leapFrog(&atoms, params.dt, box_size, cutoff, params.density);

        if (i % params.data_export_interval == 0) {
            energies[0].push_back(e.kinetic / atoms.size());
            energies[1].push_back(e.potential / atoms.size());
            energies[2].push_back(e.pressure / (atoms.size() * 3.0));

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
    energiesCols.push_back("pressure");

    exportMultiVector("energies", energiesCols, energies, 15);
    exportVector("temperature", temperatures, 15);
}

int main(int argc, char *argv[]) {

    std::cout.precision(15);
    std::cout << "Reading parameters..." << std::endl;
    Parameters params = readParameters(argv[1]);

    double width = 0.0;
    std::vector<Atom> atoms;

    if (argc > 2) {
        std::cout << "Reading atoms..." << std::endl;
        atoms = readAtomsData(params, argv[2], &width);
    } else {
        std::cout << "Initializing atoms..." << std::endl;
        atoms = initializeAtoms(params, &width);
    }

    double cutoff = pow(2.0, 1.0 / 6.0);
    std::cout << "Cut-off " << cutoff  << std::endl;
    atoms = initializeVelocities(atoms, params.temperature);


    std::cout << "Initialized " << atoms.size() << " atoms" << std::endl;

    std::cout << "Running simulation" << std::endl;
    runSimulation(params, atoms, width, cutoff);
    
    return 0;
}



#include <cmath>
#include <iostream>
#include <utility>

#include "Atoms.h"
#include "Molecules.h"
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
StepResult getPotentialAndUpdateForEach(Molecule m1, Molecule m2, double box_size, double cutoff, double b) {

    double cutoff_squared = pow(cutoff, 2.0);
    double result_potential = 0;
    double force = 0;

    Vector diff = getDistance(m1.position, m2.position, box_size, cutoff);
    double length_squared = squaredLength(diff);

    std::cout.precision(15);

    if (length_squared < cutoff_squared) {

        for (int i = 0; i < m1.sites.size(); i++) {
            for (int j = 0; j < m2.sites.size(); j++) {
                Site s1 = m1.sites[i];
                Site s2 = m2.sites[j];
                int typeSum = s1.type + s2.type;

                if (s1.type == s2.type || typeSum == 5) {

                    Vector site_diff = getDistance(s1.position, s2.position, box_size, cutoff);
                    double site_length_squared_inv = 1.0 / squaredLength(site_diff);

                    double site_potential = 0;
                    double site_force = 0;

                    switch (typeSum) {
                        case 2:
                            double site_len_inv_cubed = std::pow(site_length_squared_inv, 3.0);
                            site_potential = 4.0 * site_len_inv_cubed * (site_len_inv_cubed - 1.0);
                            site_force = (48.0 * site_len_inv_cubed * (site_len_inv_cubed - 0.5)) * site_length_squared_inv;
                            break;
                        case 4:
                            site_potential = 4.0 * b * std::sqrt(site_length_squared_inv);
                            site_force = site_potential * site_length_squared_inv;
                            break;
                        case 5:
                            site_potential = -2.0 * b * std::sqrt(site_length_squared_inv);
                            site_force = site_potential * site_length_squared_inv;
                            break;
                        case 6:
                            site_potential = b * std::sqrt(site_length_squared_inv);
                            site_force = site_potential * site_length_squared_inv;
                            break;
                    }

                    Vector pairPotential = scale(site_diff, site_force);
                    s1.force = sum(s1.force, pairPotential);
                    s2.force = subtract(s2.force, pairPotential);
                    result_potential += site_potential;
                }
            }
        }
    }

    return {result_potential, force * length_squared, m1, m2};
}

LeapFrogResult getTotalPotentialAndUpdateForEach(std::vector<Molecule> *molecules, double box_size, double cutoff) {

    double potential_energy = 0.0;
    double forces = 0.0;

    for (int i = 0; i < (*molecules).size(); i++) {
        for (int j = i + 1; j < (*molecules).size(); j++) {
            StepResult sr = getPotentialAndUpdateForEach((*molecules)[i], (*molecules)[j], box_size, cutoff);
            potential_energy += sr.potentialEnergy;
            forces += sr.force;
            (*molecules)[i] = sr.m1;
            (*molecules)[j] = sr.m2;
        }
    }

    return {potential_energy, forces};
}

Energies leapFrog(std::vector<Molecule> *molecules, double dt, double box_size, double cutoff, double density) {

    double n = (*molecules).size();
    double half_dt = dt / 2.0;

    for (int i = 0; i < n; i++) {
        Molecule m = (*molecules)[i];
        m.velocity = sum(m.velocity, scale(m.acceleration, half_dt));
        m.position = applyBoundaries(sum(m.position, scale(m.velocity, dt)), box_size);
        std::pair<double, double> min_max= getMinMaxCoordinate(m.position);

        if (min_max.first < 0 || min_max.second > box_size) {
            Vector p = m.position;
            std::cout << "\nMolecule " << i << " out of bounds (" << p.x << ',' << p.y << ',' << p.z << "), box size: " << box_size << std::endl;
            throw 100;
        }

        // Resetting forces
        std::vector<Site> zeroedSites;
        for (Site site : m.sites) {
            site.force = getZeroVector();
            zeroedSites.push_back(site);
        }
        m.sites = zeroedSites;

        (*molecules)[i] = m;
    }

    // calculate velocity and position for half step
    LeapFrogResult lfResult = getTotalPotentialAndUpdateForEach(molecules, box_size, cutoff);

    double kinetic_energy = 0;
    for (int i = 0; i < n; i++) {
        Molecule m = (*molecules)[i];
        m.velocity = sum(m.velocity, scale(m.acceleration, half_dt));
        (*molecules)[i] = m;
        Vector w = computeAngularVelocities(m);
        kinetic_energy += dot(inertia, dot(w, w));
    }

    double pressure = density * (kinetic_energy + lfResult.forces) / ((*molecules).size() * 3.0);

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

    std::cout << "Initializing molecules..." << std::endl;
    atoms = initializeAtoms(params, &width);

    double cutoff = std::pow(2.0, 1.0 / 6.0);
    std::cout << "Cut-off " << cutoff  << std::endl;
    atoms = initializeVelocities(atoms, params.temperature);




    std::cout << "Initialized " << atoms.size() << " atoms" << std::endl;

    std::cout << "Running simulation" << std::endl;
    runSimulation(params, atoms, width, cutoff);
    
    return 0;
}



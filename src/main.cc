#include <cmath>
#include <iostream>
#include <utility>

#include "Atoms.h"
#include "Molecules.h"
#include "Energies.h"
#include "Parameters.h"
#include "StepResult.h"
#include "PredictorCorrector.h"
#include "utils/OutputUtils.h"

// eps = 0.155 kcal/mol (1.08e-14 erg/molecule)
// theta (o) = 3.154 A
// b = e^2/eps*theta
const double B = 183.5;

struct StepEnergies {
    double totalPotential;
    double forces;
};

double applyBoundary(double position, double box_size) {
    if (position < 0) {
        return position + box_size;
    } else if (position >= box_size) {
        return position - box_size;
    }

    return position;
}

void checkOutOfBounds(Molecule m, double box_size) {
    std::pair<double, double> min_max= getMinMaxCoordinate(m.position);
    if (min_max.first < 0 || min_max.second > box_size) {
        Vector p = m.position;
        std::cout << "\nMolecule out of bounds (" << p.x << ',' << p.y << ',' << p.z << "), box size: " << box_size << std::endl;
        throw 100;
    }
}

Vector applyBoundaries(std::vector<Molecule> *molecules, double box_size) {

    for (int i = 0; i < molecules->size(); i++) {

        Molecule m = molecules->[i];

        m.position = {
                applyBoundary(m.position.x, box_size),
                applyBoundary(m.position.y, box_size),
                applyBoundary(m.position.z, box_size)
        };

        checkOutOfBounds(m, box_size);

        molecules->[i] = m;
    }
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
StepResult getPotentialAndUpdateForEach(Molecule *m1, Molecule *m2, double box_size, double cutoff, double b) {

    double cutoff_squared = pow(cutoff, 2.0);
    double result_potential = 0;
    double force = 0;

    Vector diff = getDistance(m1->position, m2->position, box_size, cutoff);
    double length_squared = squaredLength(diff);

    std::cout.precision(15);

    if (length_squared < cutoff_squared) {

        for (int i = 0; i < m1->sites.size(); i++) {
            for (int j = 0; j < m2->sites.size(); j++) {
                Site s1 = m1->sites[i];
                Site s2 = m2->sites[j];
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
                m1->sites[i] = s1;
                m2->sites[j] = s2;
            }
        }
    }

    return {result_potential, force * length_squared};
}

StepEnergies calculatePotentials(std::vector<Molecule> *molecules, double box_size, double cutoff) {

    double potential_energy = 0.0;
    double forces = 0.0;

    for (int i = 0; i < molecules->size(); i++) {
        for (int j = i + 1; j < molecules->size(); j++) {
            StepResult sr = getPotentialAndUpdateForEach(*(molecules->[i]), *(molecules->[j]), box_size, cutoff, B);
            potential_energy += sr.potentialEnergy;
            forces += sr.force;
        }
    }

    return {potential_energy, forces};
}

// Checked
void normalizeTemperature(std::vector<Molecule> *molecules) {

    double n = molecules->size();
    double top = 0;
    double bottom = 0;

    for (Molecule m : *molecules) {
        Vector w = computeAngularVelocities(m);
        top += dot(m.velocity, m.acceleration) + dot(w, m.torque);
        bottom += squaredLength(m.velocity) + dot(m.inertia, w);
    }

    double vDiff = -top/bottom;

    for (int i = 0; i < n; i++) {
        Molecule m = molecules->[i];
        m.acceleration = sum(m.acceleration, scale(m.velocity, vDiff));
        m.qAcceleration = sum(m.qAcceleration, scale(m.qVelocity, vDiff));
        molecules->[i] = m;
    }
}

void adjustTemperature(std::vector<Molecule> *molecules, double velocityScale) {
    double velocityDiff = 0;
    long n = molecules->size();
    double velocitiesSum = 0;

    for (Molecule m : *molecules) {
        Vector w = computeAngularVelocities(m);
        velocitiesSum += dot(m.inertia, w, w);
    }

    velocityDiff = velocityScale / sqrt(velocitiesSum / (double) n);

    for (int i = 0; i < n; i++) {
        Molecule m = molecules->[i];
        m.qVelocity = scale(m.qVelocity, velocityDiff);
        molecules->[i] = m;
    }
}

void adjustEquilibrationTemperature(std::vector<Molecule> *molecules, long equilibrationInterval, double velocityScale, double accKineticEnergy) {

    long n = molecules->size();
    double velocityDiff = velocityScale / sqrt(2.0 * accKineticEnergy / equilibrationInterval);

    for (int i = 0; i < n; i++) {
        Molecule m = molecules->[i];
        m.velocity = scale(m.velocity, velocityDiff);
        molecules->[i] = m;
    }
}

Energies predictorCorrectorStep(
        std::vector<Molecule> *molecules,
        double dt,
        double box_size,
        double cutoff,
        double density
) {

    predict(molecules, dt);
    predictQuaternion(molecules, dt);

    updateSitesCoordinates(molecules);

    StepEnergies stepEnergies = calculatePotentials(molecules, box_size, cutoff);
    calculateTorques(molecules);
    computeAccelerationQuats(molecules);

    normalizeTemperature(molecules);

    correct(molecules, dt);
    correctQuaternion(molecules, dt);

    adjustQuaternions(molecules);
    applyBoundaries(molecules, box_size);


    // Calculate measurements
    double kinetic_energy = 0;
    for (Molecule m : *molecules) {
        Vector w = computeAngularVelocities(m);
        kinetic_energy += dot(m.inertia, w, w);
    }

    double pressure = density * (kinetic_energy + stepEnergies.forces) / (molecules->size() * 3.0);

    return Energies(
            0.5 * kinetic_energy / molecules->size(),
            stepEnergies.totalPotential / molecules->size(),
            pressure
    );
}

void runSimulation(Parameters params, std::vector<Molecule> *molecules, double box_size, double cutoff, double velocityScale) {

    double time = 0.0;
    long i = 0;
    double kineticEnergySum = 0;

    std::vector<std::vector<double>> energies;
    std::vector<double> kinetic;
    std::vector<double> totalEnergy;
    std::vector<double> pressure;
    energies.push_back(kinetic);
    energies.push_back(totalEnergy);
    energies.push_back(pressure);


    while (time < params.max_time) {

        Energies e = predictorCorrectorStep(molecules, params.dt, box_size, cutoff, params.density);

        // Temperature adjustments
        if (i > params.equilibration_steps) {
            if (i % params.adjust_temperature_interval == 0) {
                adjustTemperature(molecules, velocityScale);
            }
        } else {
            kineticEnergySum += e.kinetic;

            if (i % params.adjust_equilibration_temp_interval == 0) {
                adjustEquilibrationTemperature(molecules, params.adjust_equilibration_temp_interval, velocityScale, kineticEnergySum);
                kineticEnergySum = 0;
            }
        }

        if (i % params.data_export_interval == 0) {
            energies[0].push_back(e.kinetic);
            energies[1].push_back(e.kinetic + e.potential);
            energies[2].push_back(e.pressure);
            std::cout << ".";
        }

        time += params.dt;
        i++;
    }

    std::cout << std::endl;

    std::cout << "Exporting results" << std::endl;
    std::vector<std::string> energiesCols;
    energiesCols.push_back("kinetic");
    energiesCols.push_back("totalEnergy");
    energiesCols.push_back("pressure");

    exportMultiVector("energies", energiesCols, energies, 15);
}

int main(int argc, char *argv[]) {

    std::cout.precision(15);
    std::cout << "Reading parameters..." << std::endl;
    Parameters params = readParameters(argv[1]);

    std::cout << "Initializing molecules..." << std::endl;

    double width = 0.0;
    double velocityScale = 0.0;

    std::vector<Molecule> molecules = initializeMolecules(
            params,
            getWaterMoleculeStub(),
            &velocityScale,
            &width
    );

    std::cout << "Initialized " << molecules.size() << " molecules" << std::endl;

    std::cout << "Running simulation" << std::endl;
    runSimulation(params, &molecules, width, params.cut_off, velocityScale);
    
    return 0;
}



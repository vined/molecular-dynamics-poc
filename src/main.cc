#include <cmath>
#include <iostream>
#include <GL/glut.h>

#include "Atoms.h"
#include "Parameters.h"


// Argon parameters
#define ARGON_SIGMA 3.4e-8 // ergs
#define ARGON_MASS 6.6e-23
#define ARGON_EPSILON 1.66e-14

// Potential energy cutoff
#define CUTOFF 2.0


// TODO:
// - initial state:
//   - initial velocities depend on temperature
//   - assign velocities according to Boltzmann or uniform distribution
//   - total system momentum should be zero
// - solvent


// Forces and Energies

double getKineticEnergy() {
    // sum(p[i]^2 / 2*m)
    // is it the same as sum(a[i]^2 / 2*m) ?
    return 0;
}

// Total potential energy
// U = sum[i=1..N]( sum[i=j>i..N]( u( |r_i - r_j| ))
double getPotentialEnergy() {
    return 0;
}

struct Vector getForce(double mass, struct Vector p1, struct Vector p2, double dt) {
    struct Vector f = {};
    return f;
}


//Lenard-Jones Potential for van der Waals system
double getPotentialAndUpdateForEach(struct Atom a1, struct Atom a2) {

    // Todo move these to parent function
    double sigma_squared = pow(ARGON_SIGMA, 2.0),
            cutoff_squared = pow(CUTOFF, 2.0),
            repulsion_erg = 48.0 * ARGON_EPSILON,
            system_repulsion_erg = 4.0 * ARGON_EPSILON;


    struct Vector diff = subtract(a1.position, a2.position);
    double length_squared = std::pow(diff.x, 2.0) + std::pow(diff.y, 2.0) + std::pow(diff.z, 2.0);

    if (length_squared < cutoff_squared) {

        // Temporary variables to speed up calculations
        double a_squared = sigma_squared / length_squared, // (sigma / r_ij) ^ 2
                a_sixth = std::pow(a_squared, 3.0); // (sigma / r_ij) ^ 2

        double potential = repulsion_erg * a_sixth * (a_sixth - 0.5);
        struct Vector pairPotential = scale(diff, potential);

        a1.potential = sum(a1.potential, pairPotential);
        a2.potential = subtract(a2.potential, pairPotential);

        // return system potential
        return system_repulsion_erg * a_sixth * (a_sixth - 1.0);
    }
    return 0;
}

double getTotalPotentialAndUpdateForEach(std::vector<Atom> atoms) {

    double potential_energy = 0.0;

    for (int i = 0; i < atoms.size(); i++) {
        for (int j = i + 1; j < atoms.size(); j++) {
            potential_energy += getPotentialAndUpdateForEach(atoms[i], atoms[j]);
        }
    }

    return potential_energy
}

// Velocity-Verlet method
std::pair velocityVerlet(std::vector<Atom> atoms, double dt) {

    double half_dt = dt / 2.0,
            mass_inv = 1 / ARGON_MASS;

    // calculate velocity and position for half step
    for (Atom a : atoms) {
        a.velocity = sum(a.velocity, scale(a.acceleration, half_dt));
        a.position = sum(a.position, scale(a.velocity, dt));
    }

    double potential_energy = getTotalPotentialAndUpdateForEach(atoms),
            kinetic_energy = 0;

    for (Atom a : atoms) {
        a.acceleration = scale(a.potential, mass_inv);
        a.velocity = sum(a.velocity, scale(a.acceleration, half_dt));
        kinetic_energy += squaredLength(a.velocity);
    }

    kinetic_energy = ARGON_MASS * kinetic_energy / 2.0;
    return {kinetic_energy, potential_energy};
}

void runSimulation(Parameters parameters) {

    double time = 0;
    std::vector<Atom> atoms = parameters.atoms;

    while (time < parameters.max_time) {
        std::pair energies = velocityVerlet(atoms, parameters.dt);
    }
}

void keyboardFn(unsigned char key, int x, int y) {
    switch (key) {
        case 'p':
        case ' ':
//            pause();
            glutPostRedisplay();
            break;
        case 'x':
            exit(0);
    }
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}



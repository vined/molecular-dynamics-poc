#include <cmath>
#include <iostream>

#include "vector-utils.h"
#include "parameters.h"

// Argon parameters
#define ARGON_SIGMA 3.4e-8 // ergs
#define ARGON_MASS 6.6e-23
#define ARGON_EPSILON 1.66e-14

// Potential energy cutoff
#define CUTOFF 2.0

struct AtomType {
    int type;
    double mass;
};

struct Atom {
    struct AtomType type;
    struct Vector position;
    struct Vector potential;
    struct Vector acceleration;
};


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
double getTotalPotentialAndUpdateForEach(struct Atom a1, struct Atom a2) {
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

        double inter_potential = repulsion_erg * a_sixth * (a_sixth - 0.5);
        a1.potential = add(a1.potential, inter_potential);
        a2.potential = add(a2.potential, -inter_potential);

        double rij = length(subtract(a1.position, a2.position)),
                div = ARGON_SIGMA / rij,
                repulsion = pow(div, 12.0),
                attraction = pow(div, 6);

        // update atom potentials


        // return system potential
        return system_repulsion_erg * a_sixth * (a_sixth - 0.5);
    }
    return 0;
}


// Velocity-Verlet method
// 1. get next acceleration
// 2. get next velocity
// 3. get next position

// Newtown equation F = m * a
struct Vector getNextAcceleration(double mass, struct Vector force) {
    return scale(force, 1.0 / mass);
}

// v(t + dt) = v(t) + dt * (a(t) + a(t + dt)) / 2
struct Vector getNextVelocity(
        struct Vector velocity,
        struct Vector accelerationPrev,
        struct Vector accelerationNext,
        double dt
) {
    return sum(
            velocity,
            scale(sum(accelerationPrev, accelerationNext), dt / 2)
    );
}

// r(t + dt) = r(t) + dt * v(t) + dt^2 * a(t) / 2
struct Vector getNextPosition(
        struct Vector position,
        struct Vector velocity,
        struct Vector acceleration,
        double dt
) {
    return sum(
            position,
            scale(velocity, dt),
            scale(acceleration, pow(dt, 2.0) / 2)
    );
}


// Trajectories

double getAcceleration(double dv, double dt) {
    return dv / dt;
}

double getVelocity(double dx, double dt) {
    return dx / dt;
}

double getTraveledDistance(double a, double dt, double v0, double x0) {
    return a * std::pow(dt, 2.0) / 2.0 + v0 * dt + x0;
}

int main() {
    std::cout << "Hello, World!" << std::endl;

    return 0;
}



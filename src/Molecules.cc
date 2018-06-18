#include <cmath>
#include <random>
#include <iostream>

#include "Molecules.h"
#include "utils/OutputUtils.h"


// Checked
Quaternion eulerToQuaternion(std::vector<double> eAngular) {

    double a1 = 0.5 * eAngular[1];
    double a2 = 0.5 * (eAngular[0] - eAngular[2]);
    double a3 = 0.5 * (eAngular[0] + eAngular[2]);
    return {
            std::sin(a1) * std::cos(a2),
            std::sin(a1) * std::sin(a2),
            std::cos(a1) * std::cos(a3),
            std::cos(a1) * std::cos(a3),
    };
}

Molecule getWaterMoleculeStub() {
    Molecule m = Molecule();

    m.sites.push_back({1, getZeroVector(), {0, 0, -0.0206}}); // O
    m.sites.push_back({2, getZeroVector(), {0, 0, 0.0274}}); // mass center
    m.sites.push_back({3, getZeroVector(), {0, 0.24, 0.165}}); // H
    m.sites.push_back({3, getZeroVector(), {0, -0.24, 0.165}}); // H

    m.inertia = {0.0098, 0.0034, 0.0064};

    return m;
}

// Checked
void initializeAngularCoordinates(std::vector<Molecule> *molecules) {

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::default_random_engine re;


    for (int i = 0; i < molecules->size(); i++) {

        Molecule m = molecules->at(i);
        Vector random = {unif(re), unif(re), unif(re)};

        std::vector<double> eAngular;
        eAngular.push_back(std::atan2(random.x, random.y));
        eAngular.push_back(std::acos(random.z));
        eAngular.push_back(2.0 * M_PI * unif(re));

        m.quaternion = eulerToQuaternion(eAngular);

        molecules->at(i) = m;
    }
}

void initializeAngularVelocities(std::vector<Molecule> *molecules, double velocityScale) {

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::default_random_engine re;

    for (int i = 0; i < molecules->size(); i++) {

        Molecule m = molecules->at(i);
        Vector r = {unif(re), unif(re), unif(re)};

        Quaternion q = {r.x, r.y, r.z, 0};
        double f = 0.5 * velocityScale / sqrt(dot(m.inertia, r, r));

        m.qVelocity = scale(multiply(m.quaternion, q), f);

        molecules->at(i) = m;
    }
}

void initializeVelocities(std::vector<Molecule> *molecules, double velocityScale) {

    long n = molecules->size();
    Vector totalVelocity = getZeroVector();

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::default_random_engine re;

    for (int i = 0; i < n; i++) {
        Molecule m = molecules->at(i);

        Vector r = {unif(re), unif(re), unif(re)};
        m.velocity = scale(r, velocityScale);

        molecules->at(i) = m;
        totalVelocity = sum(totalVelocity, m.velocity);
    }

    Vector velocityReducer = scale(totalVelocity, -1.0 / (double) n);

    for (int i = 0; i < n; i++) {
        Molecule m = molecules->at(i);
        m.velocity = sum(m.velocity, velocityReducer);
        molecules->at(i) = m;
    }
}

void resetStepValues(std::vector<Molecule> *molecules) {
    long n = molecules->size();

    for (int i = 0; i < n; i++) {
        Molecule m = molecules->at(i);

        m.hydrogenBonds = 0;

        molecules->at(i) = m;
    }
}

// Checked
Vector computeAngularVelocities(Molecule m) {

    Quaternion q = m.qVelocity;
    q.d = -q.d;

    Quaternion result = scale(multiply(q, m.quaternion), 2.0);

    return {
            result.a,
            result.b,
            result.c
    };
}

std::vector<Molecule> generateTwoMolecules(Parameters params, Molecule stub, double *width) {

    std::vector<Molecule> molecules;

    *width = params.density * 2.0;
    double half_step = params.density / 2.0;

    Molecule m1 = stub;
    m1.position = {half_step, half_step, half_step};
    molecules.push_back(m1);

    Molecule m2 = stub;
    m2.position = {*width - half_step, *width - half_step, *width - half_step};
    molecules.push_back(m2);

    return molecules;
}

std::vector<Molecule> generateMoleculesWithPositions(Parameters params, Molecule stub, double *width) {

    std::vector<Molecule> molecules;

    double half_step = params.density / 2.0;
    double half_half_step = half_step / 2.0;
    double x = half_half_step;

    for (int i = 0; i < params.lattice_count; i++) {

        double y = 0;
        for (int j = 0; j < params.lattice_count; j++) {

            double z = 0;
            for (int k = 0; k < params.lattice_count; k++) {

                Molecule m = stub;
                m.position = {x + half_step, y + half_step, z + half_step};
                molecules.push_back(m);

                z += params.density;
            }
            y += params.density;
        }
        x += params.density;
    }

    *width = x - half_half_step;

    return molecules;
}

std::vector<Molecule> initializeMolecules(Parameters params, Molecule stub, double *velocityScale, double *width) {

    std::vector<Molecule> molecules = generateTwoMolecules(params, stub, width);
//    std::vector<Molecule> molecules = generateMoleculesWithPositions(params, stub, width);

    *velocityScale = std::sqrt(3.0 * (1.0 - 1.0 / (double) molecules.size()) * params.temperature);
//    std::cout << "Temperature: -- " << params.temperature << " -- " << std::endl;
//    std::cout << "Before sqrt: -- " << 3.0 * (1.0 - 1.0 / (double) molecules.size()) * params.temperature << " -- " << std::endl;
//    std::cout << "Velocity scale: -- " << *velocityScale << " -- " << std::endl;

    initializeVelocities(&molecules, *velocityScale);
    initializeAngularCoordinates(&molecules);
    initializeAngularVelocities(&molecules, *velocityScale);

//    long i = 0;
//    for (Molecule m : molecules) {
//        std::cout << "Molecule in init: -- " << i << " -- " << std::endl;
//        std::cout << "Molecule position: " << vectorToString(m.position) << std::endl;
//        std::cout << "Molecule velocity: " << vectorToString(m.velocity) << std::endl;
//        std::cout << "Molecule inertia: " << vectorToString(m.inertia) << std::endl;
//        std::cout << "Molecule torque: " << vectorToString(m.torque) << std::endl;
//        std::cout << "Molecule angular coords: " << quatToString(m.quaternion) << std::endl;
//        std::cout << "Molecule angular velocities: " << quatToString(m.qVelocity) << std::endl;
//        i++;
//    }

    return molecules;
}

// Checked
void computeAccelerationQuats(std::vector<Molecule> *molecules) {

    for (int i = 0; i < molecules->size(); i++) {

        Molecule m = molecules->at(i);
        Vector w = computeAngularVelocities(m);

//        std::cout << "W: " << vectorToString(w) << std::endl;
//        std::cout << "Inertia: " << vectorToString(m.inertia) << std::endl;

        Quaternion q = {
                (m.torque.x + (m.inertia.y - m.inertia.z) * w.y * w.z) / m.inertia.x,
                (m.torque.y + (m.inertia.z - m.inertia.x) * w.z * w.x) / m.inertia.y,
                (m.torque.z + (m.inertia.x - m.inertia.y) * w.x * w.y) / m.inertia.z,
                -2.0 * squareLength(m.qVelocity)
        };

//        std::cout << "Q: " << quatToString(q) << std::endl;

        m.qAcceleration = scale(multiply(m.quaternion, q), 0.5);

//        std::cout << "Acc: " << quatToString(m.qAcceleration) << std::endl;

        molecules->at(i) = m;
    }
}

// Checked
Matrix makeRotationMatrix(Quaternion q, int transpose) {

    std::vector<double> torque = toVector(q);
    std::vector<double> p;

    for (int i = 0; i < 4; i ++) {
        for (int j = i; j < 4; j ++) {
            p.push_back(2.0 * torque[i] * torque[j]);
        }
    }

    std::vector<double> matrix;
    double sign = transpose ? 1.0 : -1.0;
    
    matrix.push_back(p[0] + p[9] - 1.0);
    matrix.push_back(p[1] + sign * p[8]);
    matrix.push_back(p[2] - sign * p[6]);
    matrix.push_back(p[1] - sign * p[8]);
    matrix.push_back(p[4] + p[9] - 1.0);
    matrix.push_back(p[5] + sign * p[3]);
    matrix.push_back(p[2] + sign * p[6]);
    matrix.push_back(p[5] - sign * p[3]);
    matrix.push_back(p[7] + p[9] - 1.0);

    return {matrix};
}

// Checked (is mol.ra == m.acceleration?)
void calculateTorques(std::vector<Molecule> *molecules) {

    for (int i = 0; i < molecules->size(); i++) {

        Molecule m = molecules->at(i);
        Vector totalTorque = getZeroVector();

        m.acceleration = getZeroVector();

        for (Site site : m.sites) {
            m.acceleration = sum(m.acceleration, site.force);

            Vector diff = subtract(site.position, m.position);
            Vector torque = crossProduct(diff, site.force);
            totalTorque = sum(totalTorque, torque);
        }

        Matrix rotationMatrix = makeRotationMatrix(m.quaternion, 0);
        m.torque = multiplyVector(rotationMatrix, totalTorque);

        molecules->at(i) = m;
    }
}

// Checked
// Initializing sites for each molecule from its position and site stub
void createSites(std::vector<Molecule> *molecules, std::vector<Site> stubSites) {

    for (int i = 0; i < molecules->size(); i++) {

        Molecule m = molecules->at(i);
        Matrix rotation = makeRotationMatrix(m.quaternion, 1);

        // Just in case it is used on mol with sites
        m.sites.clear();

        for (Site site : stubSites) {
            Vector v = multiplyVector(rotation, site.position);
            m.sites.push_back({site.type, getZeroVector(), sum(m.position, v)});
        }

        molecules->at(i) = m;
    }
}

// Checked
void updateSitesCoordinates(std::vector<Molecule> *molecules) {

    for (int i = 0; i < molecules->size(); i++) {

        Molecule m = molecules->at(i);
        Matrix rotation = makeRotationMatrix(m.quaternion, 1);

        for (int j = 0; j < m.sites.size(); j++) {
            Site site = m.sites[j];
            Vector v = multiplyVector(rotation, site.position);
            site.position = sum(m.position, v);
            m.sites[j] = site;
        }

        molecules->at(i) = m;
    }
}

// Checked
// Preventing accumulation of computation errors
void adjustQuaternions(std::vector<Molecule> *molecules) {

    for (int i = 0; i < molecules->size(); i++) {

        Molecule m = molecules->at(i);

        m.quaternion = scale(m.quaternion, 1.0 / squareLength(m.quaternion));

        molecules->at(i) = m;
    }
}

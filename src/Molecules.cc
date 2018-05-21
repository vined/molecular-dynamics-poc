#include <cmath>
#include <random>

#include "Molecules.h"


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

std::vector<Molecule> initializeAngularCoordinates(std::vector<Molecule> *molecules) {

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::default_random_engine re;


    for (int i = 0; i < (*molecules).size(); i++) {

        Molecule m = (*molecules)[i];
        Vector random = {unif(re), unif(re), unif(re)};

        std::vector<double> eAngular;
        eAngular.push_back(std::atan2(random.x, random.y));
        eAngular.push_back((double) std::acos(random.z));
        eAngular.push_back(2.0 * M_PI * unif(re));

        m.quaternion = eulerToQuaternion(eAngular);

        (*molecules)[i] = m;
    }
}

std::vector<Molecule> initializeAngularVelocities(std::vector<Molecule> *molecules, Vector inertia, double temperature) {

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::default_random_engine re;

    double velocityScale = std::sqrt(3.0 * (1.0 - (1.0 / (*molecules).size()) * temperature));

    for (int i = 0; i < (*molecules).size(); i++) {

        Molecule m = (*molecules)[i];
        Vector r = {unif(re), unif(re), unif(re)};

        Quaternion q = {r.x, r.y, r.z, 0};
        double f = 0.5 * velocityScale / length(dot(inertia, dot(r, r)));

        m.qVelocity = scale(multiply(m.quaternion, q), f);

        (*molecules)[i] = m;
    }
}

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

void computeAccelerationQuats(std::vector<Molecule> *molecules, Vector inertia) {

    for (int i = 0; i < (*molecules).size(); i++) {

        Molecule m = (*molecules)[i];
        Vector w = computeAngularVelocities(m);

        Quaternion q = {
                (m.torque.x + (inertia.y - inertia.z) * w.y * w.z) / inertia.x,
                (m.torque.y + (inertia.z - inertia.x) * w.z * w.x) / inertia.y,
                (m.torque.z + (inertia.x - inertia.y) * w.x * w.y) / inertia.z
        };

        m.qAcceleration = scale(multiply(m.quaternion, q), 0.5);

        (*molecules)[i] = m;
    }
}

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

void calculateTorques(std::vector<Molecule> *molecules) {

    for (int i = 0; i < (*molecules).size(); i++) {

        Molecule m = (*molecules)[i];
        Vector totalTorque = getZeroVector();

        for (Site site : m.sites) {
            m.position = sum(m.position, site.force);

            Vector diff = subtract(site.position, m.position);
            Vector torque = crossProduct(diff, site.force);
            totalTorque = sum(totalTorque, torque);
        }

        Matrix rotationMatrix = makeRotationMatrix(m.quaternion, 0);
        m.torque = multiplyVector(rotationMatrix, totalTorque);

        (*molecules)[i] = m;
    }
}

void initializeSites(std::vector<Molecule> *molecules, std::vector<Site> sites) {

    for (int i = 0; i < (*molecules).size(); i++) {

        Molecule m = (*molecules)[i];
        Matrix rotation = makeRotationMatrix(m.quaternion, 1);

        for (Site site : sites) {

            Vector v = multiplyVector(rotation, site.position);

            m.sites.push_back({site.type, getZeroVector(), sum(m.position, v)});
        }

        (*molecules)[i] = m;
    }
}

void adjustQuaternions(std::vector<Molecule> *molecules) {

    for (int i = 0; i < (*molecules).size(); i++) {

        Molecule m = (*molecules)[i];

        m.quaternion = scale(m.quaternion, length(m.quaternion));

        (*molecules)[i] = m;
    }
}

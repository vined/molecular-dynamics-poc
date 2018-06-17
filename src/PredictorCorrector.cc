

#include "PredictorCorrector.h"
#include "Vectors.h"


const double divisor = 24.0;
const double predictorPositionCoefs[] = {19.0, -10.0, 3.0};
const double predictorVelocityCoefs[] = {27.0, -22.0, 7.0};
const double correctorPositionCoefs[] = {3.0, 10.0, -1.0};
const double correctorVelocityCoefs[] = {7.0, 6.0, -1.0};


Vector fromPreviousStep(Vector acceleration, Vector acceleration1, Vector acceleration2, double w, double const coefs[]) {

    Vector v1 = scale(acceleration, coefs[0]);
    Vector v2 = scale(acceleration1, coefs[1]);
    Vector v3 = scale(acceleration2, coefs[2]);

    return scale(sum(v1, sum(v2, v3)), w);

}

Quaternion fromPreviousStepQuat(Quaternion qAcceleration, Quaternion qAcceleration1, Quaternion qAcceleration2, double w, double const coefs[]) {

    Quaternion q1 = scale(qAcceleration, coefs[0]);
    Quaternion q2 = scale(qAcceleration1, coefs[1]);
    Quaternion q3 = scale(qAcceleration2, coefs[2]);

    return scale(sum(q1, sum(q2, q3)), w);

}

Vector pcPosition4(
        Vector positionPrev, Vector velocity,
        Vector acceleration, Vector acceleration1, Vector acceleration2,
        double dt, double wPosition, double const positionCoefs[]) {

    Vector v1 = scale(velocity, dt);
    Vector v2 = fromPreviousStep(acceleration, acceleration1, acceleration2, wPosition, positionCoefs);

    return sum(positionPrev, sum(v1, v2));
}

Quaternion pcPosition4Quaternion(
        Quaternion quaternionPrev, Quaternion qVelocity,
        Quaternion qAcceleration, Quaternion qAcceleration1, Quaternion qAcceleration2,
        double dt, double wPosition, double const positionCoefs[]) {

    Quaternion q1 = scale(qVelocity, dt);
    Quaternion q2 = fromPreviousStepQuat(qAcceleration, qAcceleration1, qAcceleration2, wPosition, positionCoefs);

    return sum(quaternionPrev, sum(q1, q2));
}

Vector pcVelocity4(
        Vector position, Vector positionPrev,
        Vector acceleration, Vector acceleration1, Vector acceleration2,
        double dt, double wVelocity, double const velocityCoefs[]) {

    Vector v1 = scale(subtract(position, positionPrev), 1.0 / dt);
    Vector v2 = fromPreviousStep(acceleration, acceleration1, acceleration2, wVelocity, velocityCoefs);

    return sum(v1, v2);
}

Quaternion pcVelocity4Quaternion(
        Quaternion quaternion, Quaternion quaternionPrev,
        Quaternion qAcceleration, Quaternion qAcceleration1, Quaternion qAcceleration2,
        double dt, double wVelocity, double const velocityCoefs[]) {

    Quaternion q1 = scale(subtract(quaternion, quaternionPrev), 1.0 / dt);
    Quaternion q2 = fromPreviousStepQuat(qAcceleration, qAcceleration1, qAcceleration2, wVelocity, velocityCoefs);

    return sum(q1, q2);
}

void predict(std::vector<Molecule> *molecules, double dt) {

    double wPosition = pow(dt, 2.0) / divisor;
    double wVelocity = dt / divisor;

    for (int i = 0; i < molecules->size(); i++) {

        Molecule m = molecules->[i];

        m.positionPrev = m.position;
        m.velocityPrev = m.velocity;

        m.position = pcPosition4(m.position, m.velocity, m.acceleration, m.acceleration1, m.acceleration2, dt, wPosition, predictorPositionCoefs);
        m.velocity = pcVelocity4(m.position, m.positionPrev, m.acceleration, m.acceleration1, m.acceleration2, dt, wVelocity, predictorVelocityCoefs);

        m.acceleration2 = m.acceleration1;
        m.acceleration1 = m.acceleration;

        molecules->[i] = m;
    }
}

void predictQuaternion(std::vector<Molecule> *molecules, double dt) {

    double wPosition = pow(dt, 2.0) / divisor;
    double wVelocity = dt / divisor;

    for (int i = 0; i < molecules->size(); i++) {

        Molecule m = molecules->[i];

        m.quaternionPrev = m.quaternion;
        m.qVelocityPrev = m.qVelocity;

        m.quaternion = pcPosition4Quaternion(m.quaternionPrev, m.qVelocity, m.qAcceleration, m.qAcceleration1, m.qAcceleration2, dt, wPosition, predictorPositionCoefs);
        m.qVelocity = pcVelocity4Quaternion(m.quaternion, m.quaternionPrev, m.qAcceleration, m.qAcceleration1, m.qAcceleration2, dt, wVelocity, predictorVelocityCoefs);

        m.qAcceleration2 = m.qAcceleration1;
        m.qAcceleration1 = m.qAcceleration;

        molecules->[i] = m;
    }
}

void correct(std::vector<Molecule> *molecules, double dt) {

    double wPosition = pow(dt, 2.0) / divisor;
    double wVelocity = dt / divisor;

    for (int i = 0; i < molecules->size(); i++) {

        Molecule m = molecules->[i];

        m.position = pcPosition4(m.positionPrev, m.velocityPrev, m.acceleration, m.acceleration1, m.acceleration2, dt, wPosition, correctorPositionCoefs);
        m.velocity = pcVelocity4(m.position, m.positionPrev, m.acceleration, m.acceleration1, m.acceleration2, dt, wVelocity, correctorVelocityCoefs);

        molecules->[i] = m;
    }
}

void correctQuaternion(std::vector<Molecule> *molecules, double dt) {

    double wPosition = pow(dt, 2.0) / divisor;
    double wVelocity = dt / divisor;

    for (int i = 0; i < molecules->size(); i++) {

        Molecule m = molecules->[i];

        m.quaternion = pcPosition4Quaternion(m.quaternionPrev, m.qVelocityPrev, m.qAcceleration, m.qAcceleration1, m.qAcceleration2, dt, wPosition, correctorPositionCoefs);
        m.qVelocity = pcVelocity4Quaternion(m.quaternion, m.quaternionPrev, m.qAcceleration, m.qAcceleration1, m.qAcceleration2, dt, wVelocity, correctorVelocityCoefs);

        molecules->[i] = m;
    }
}

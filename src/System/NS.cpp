//
// Created by Indra on 18/09/20.
//

#include "NS.h"

NS::NS() {
    eqs = 4;
    gamma = 1.4;
    gamma_1 = gamma - 1.0;
    R = 8.3145;
    qnc = 0.0, ql_ = 0.0, grcsq = 0.0, grcsqu = 0.0, grcsqv = 0.0,
    grcsqqsq = 0.0, qnrc = 0.0, nxrc = 0.0, nyrc = 0.0;
    leftDensitySqrt = 1.0, rightDensitySqrt = 1.0;
}

void NS::setConstant(const double &SpecificHeat, const double &GasConstant) {
    eqs = 4;
    gamma = SpecificHeat;
    gamma_1 = gamma - 1.0;
    R = GasConstant;
}

int NS::numberOfEquations() const {
    return eqs;
}

void NS::primitive(const Vector4d &q, Vector4d &pv) const {
    pv(rho) = q(rho);
    pv(vex) = q(momentx) / q(rho);
    pv(vey) = q(momenty) / q(rho);
    pv(p) = (q(e) - 0.5 * (q(momentx) * pv(vex) + q(momenty) * pv(vey))) * gamma_1;
}

void NS::conservative(const Vector4d &pv, Vector4d &q) const {
    q(rho) = pv(rho);
    q(momentx) = pv(rho) * pv(vex);
    q(momenty) = pv(rho) * pv(vey);
    q(e) = pv(p) / gamma_1 + 0.5 * (q(momentx) * pv(vex) + q(momenty) * pv(vey));
}

void NS::Enthalpy(const Vector4d &pv, const Vector4d &q, double &h) {
    h = (q(e) + pv(p)) / pv(rho);
}

void NS::RoeAveraging(const double &leftEnthalpy, const double &rightEnthalpy, const Vector4d &leftpv,
                      const Vector4d &rightpv, Vector2d &meanVelocity, double &meanEnthalpy, double &meanDynPressure,
                      double &soundSpd) {

    leftDensitySqrt = sqrt(fabs(leftpv(rho)));
    rightDensitySqrt = sqrt(fabs(rightpv(rho)));
    meanDensity(0) = leftDensitySqrt / (leftDensitySqrt + rightDensitySqrt);
    meanDensity(1) = 1.0 - meanDensity(0);
    meanVelocity = meanDensity(0) * leftpv.segment<2>(vex) + meanDensity(1) * rightpv.segment<2>(vex);
    meanEnthalpy = meanDensity(0) * leftEnthalpy + meanDensity(1) * rightEnthalpy;
    meanDynPressure = 0.5 * (meanVelocity(0) * meanVelocity(0) + meanVelocity(1) * meanVelocity(1));
    soundSpd = sqrt(gamma_1 * fabs(meanEnthalpy - meanDynPressure));
}

void NS::EigenVectorRight(const Vector2d &velocity, const double &enthalpy, const double &dynPressureSq,
                          const double &soundSpd, const Vector2d &norm, Matrix4d &RightEV) {
    soundSpdNorm = soundSpd * norm;
    qnc = velocity.dot(norm) * soundSpd;
    ql_ = -velocity(0) * norm(1) + velocity(1) * norm(0);

    RightEV.block<1, 4>(0, 0) << 1.0, 0.0, 1.0, 1.0;
    RightEV.block<1, 4>(1, 0) << velocity(0) - soundSpdNorm(0), -norm(1), velocity(0), velocity(0) + soundSpdNorm(0);
    RightEV.block<1, 4>(2, 0) << velocity(1) - soundSpdNorm(1), norm(0), velocity(1), velocity(1) + soundSpdNorm(1);
    RightEV.block<1, 4>(3, 0) << enthalpy - qnc, ql_, dynPressureSq, enthalpy + qnc;
}

void NS::EigenVectorLeft(const Vector2d &velocity, const double &dynPressureSq, const double &soundSpd,
                         const Vector2d &norm, Matrix4d &LeftEV) {
    grcsq = gamma_1 / soundSpd / soundSpd;
    grcsqu = grcsq * velocity(0);
    grcsqv = grcsq * velocity(1);
    grcsqqsq = dynPressureSq * grcsq;
    qnrc = velocity.dot(norm) / soundSpd;
    ql_ = -velocity(0) * norm(1) + velocity(1) * norm(0);
    nxrc = norm(0) / soundSpd;
    nyrc = norm(1) / soundSpd;

    LeftEV.block<1, 4>(0, 0) << 0.5 * (grcsqqsq + qnrc), -0.5 * (grcsqu + nxrc),
            -0.5 * (grcsqv + nyrc), 0.5 * grcsq;
    LeftEV.block<1, 4>(1, 0) << -ql_, -norm(1), norm(0), 0.0;
    LeftEV.block<1, 4>(2, 0) << 1.0 - grcsqqsq, grcsqu, grcsqv, -grcsq;
    LeftEV.block<1, 4>(3, 0) << 0.5 * (grcsqqsq - qnrc), -0.5 * (grcsqu - nxrc),
            -0.5 * (grcsqv - nyrc), 0.5 * grcsq;
}

Vector4d NS::Flux(const Vector4d &q, const Vector4d &pv, const int &direction) {
    if (direction == x) {
        fluxes(0) = q(rho) * pv(vex);
        fluxes(1) = q(momentx) * pv(vex) + pv(p);
        fluxes(2) = q(momentx) * pv(vey);
        fluxes(3) = pv(vex) * (q(e) + pv(p));
    } else if (direction == y) {
        fluxes(0) = q(rho) * pv(vey);
        fluxes(1) = q(momenty) * pv(vex);
        fluxes(2) = q(momenty) * pv(vey) + pv(p);
        fluxes(3) = pv(vey) * (q(e) + pv(p));
    }

    return fluxes;
}

Vector4d NS::maxLocalWave(const Vector4d &ql, const Vector4d &pvl, const Vector4d &qr,
                          const Vector4d &pvr, const int &direction) {
    velocity_(0) = pvl(direction + 1);
    velocity_(1) = pvr(direction + 1);
    soundSpd_(0) = sqrt(fabs(gamma * pvl(p) / pvl(rho)));
    soundSpd_(1) = sqrt(fabs(gamma * pvr(p) / pvr(rho)));

    waveSpd(0) = fmax(fabs(velocity_(0) - soundSpd_(0)), fabs(velocity_(1) - soundSpd_(1)));
    waveSpd(1) = velocity_.cwiseAbs().maxCoeff();
    waveSpd(2) = waveSpd(1);
    waveSpd(3) = fmax(fabs(velocity_(0) + soundSpd_(0)), fabs(velocity_(1) + soundSpd_(1)));

    return waveSpd;
}

Vector4d NS::RusanovFlux(const Vector4d &ql, const Vector4d &pvl, const Vector4d &qr,
                         const Vector4d &pvr, const Matrix4d &LeftEV, const Matrix4d &RightEV,
                         const int &direction) {
    return 0.5 * (Flux(ql, pvl, direction) + Flux(qr, pvr, direction)
                  - RightEV * maxLocalWave(ql, pvl, qr, pvr, direction).cwiseProduct(LeftEV * (qr - ql)));
//    return 0.5 * (Flux(ql, pvl, direction) + Flux(qr, pvr, direction));
}

void NS::IdealGasEOS(const Vector4d &pv, double &Temperature) const {
    Temperature = pv(p) / pv(rho) / R;
}




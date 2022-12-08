//
// Created by mf on 18/09/20.
//

#ifndef FREEFLOW_NS_H
#define FREEFLOW_NS_H

#include "../include/CoordinateSystems.h"
#include "../include/EigenType.h"
#include <iostream>
#include "VarName.h"
#include <vector>

class NS {
public:
    NS();
    void setConstant(const double &SpecificHeat, const double &GasConstant);
    [[nodiscard]] int numberOfEquations() const;
    void primitive(const Vector4d &q, Vector4d &pv) const;
    void IdealGasEOS(const Vector4d &pv, double &Temperature) const;
    void conservative(const Vector4d &pv, Vector4d &q) const;
    static void Enthalpy(const Vector4d &pv, const Vector4d &q, double &h);
    void RoeAveraging(const double &leftEnthalpy, const double &rightEnthalpy, const Vector4d &leftpv,
                      const Vector4d &rightpv, Vector2d &meanVelocity, double &meanEnthalpy,
                      double &meanDynPressure, double &soundSpd);
    void EigenVectorRight(const Vector2d &velocity, const double &enthalpy, const double &dynPressureSq,
                          const double &soundSpd, const Vector2d &norm, Matrix4d &RightEV);
    void EigenVectorLeft(const Vector2d &velocity, const double &dynPressureSq, const double &soundSpd,
                         const Vector2d &norm, Matrix4d &LeftEV);
    Vector4d RusanovFlux(const Vector4d &ql, const Vector4d &pvl, const Vector4d &qr,
                     const Vector4d &pvr, const Matrix4d &LeftEV, const Matrix4d &RightEV,
                     const int &direction);

private:
    Vector4d Flux(const Vector4d& q, const Vector4d &pv, const int &direction);
    Vector4d maxLocalWave(const Vector4d &ql, const Vector4d &pvl, const Vector4d &qr,
                          const Vector4d &pvr, const int &direction);

    int eqs;
    double gamma, gamma_1, R{};

    double qnc, ql_, grcsq, grcsqu, grcsqv, grcsqqsq, qnrc, nxrc, nyrc;
    double leftDensitySqrt, rightDensitySqrt;
    Vector2d meanDensity, soundSpdNorm, soundSpd_, velocity_;
    Vector4d fluxes, waveSpd;
};


#endif //FREEFLOW_NS_H

//mf
// Created by Indra on 17/09/20.
//

#ifndef FREEFLOW_DATA_H
#define FREEFLOW_DATA_H

#include "EigenType.h"
#include <vector>

struct faceData {
    Vector4d ql, qr;
    Matrix<double, 4, 2> flux;
};

struct cellData {
    Vector4d pv, q, residu;
    Vector2d pvt, qt;
    double T, enthalpy;
    std::vector<Matrix4d> RightEV, LeftEV;
    std::vector<faceData> face;
    std::vector<Vector4d> qn;
};

struct Data {
    double cfl, dt, t, tmax;
    double gamma, gasConstant;

    std::vector<cellData> cell;
};

#endif //FREEFLOW_DATA_H

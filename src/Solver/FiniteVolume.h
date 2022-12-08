//
// Created by mf on 17/09/20.
//

#ifndef FREEFLOW_FINITEVOLUME_H
#define FREEFLOW_FINITEVOLUME_H

#include "../include/Data.h"
#include "../Mesh/MeshHandler.h"
#include "../Schemes/Interpolator.h"
#include "../Schemes/SSPRK.h"
#include "../System/k_Omega.h"
#include "../System/NS.h"
#include <fstream>
#include <iomanip>

class FiniteVolume {
public:
    FiniteVolume(const std::string &filename, const std::string &userCase, const int &stencilWidth,
                 const int &numOfPeriodicBc, const double &tmax);
    void setConstant(const double &SpecificHeat, const double &gasConstant, const double &userCFL);
    void run();
    void save(const std::string& filename);
    ~FiniteVolume();
    void initialize();
    void boundary();
    void updateVariable();
    void interpolate();
    void residual();
    void timeStepping();

private:


    std::string cases;
    boost::circular_buffer<int> circularNumber;
    std::vector<Vector2d> naturalCoord, faceNorm, dnat;
    Matrix<double, Dynamic, Dynamic> qs;
    Matrix<double, Dynamic, 1> hs;

    Matrix<std::size_t, Dynamic, 1> stencils;

    MeshHandler::meshDefinition mesh;
    Data data;
    NS NavierStokes;

    Vector2d meanVel, globalMaxWaveSpd;
    double meanEnthalpy, meanDynPressure, soundSpd;
    Vector4d pvl, pvr;

    std::size_t neighbor;
    Matrix4d res;

    int stage;
    double waveSpd;

    Interpolator schemes;
    k_Omega kw;
    SSPRK ssprk;
};


#endif //FREEFLOW_FINITEVOLUME_H

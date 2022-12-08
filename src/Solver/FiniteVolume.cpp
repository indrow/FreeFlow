//
// Created by Indra on 17/09/20.
//

#include "FiniteVolume.h"

FiniteVolume::FiniteVolume(const std::string &filename, const std::string &userCase, const int &stencilWidth,
                           const int &numOfPeriodicBc, const double &tmax) {
    std::cout << "Finite volume class has been called.\n";

    /* Reading the input mesh */
    MeshHandler *mshObj;
    mshObj = new MeshHandler();
    mshObj->readMesh(filename, stencilWidth, numOfPeriodicBc, mesh);
    delete mshObj;

    /* Assign cases if selected */
    cases = userCase;

    /* Create circular faceNumber/node numbering */
    circularNumber.resize(mesh.numCellFace, 0);
    for (int i = 0; i < mesh.numCellFace; ++i) {
        circularNumber[i] = i;
    }

    /* Default value of the variable */
    meanEnthalpy = 0.0;
    meanDynPressure = 0.0;
    soundSpd = 0.0;
    neighbor = 0, stage = 1;
    waveSpd = 0.0;

    /* Define the natural coordinate system */
    naturalCoord.resize(mesh.numCellFace);
    dnat = naturalCoord;
    faceNorm.resize(mesh.numCellFace);

    naturalCoord[0] << -1.0, -1.0;
    naturalCoord[1] << 1.0, -1.0;
    naturalCoord[2] << 1.0, 1.0;
    naturalCoord[3] << -1.0, 1.0;

    for (int i = 0; i < mesh.numCellFace; ++i) {
        dnat[i] = naturalCoord[circularNumber[i + 1]] - naturalCoord[i];
        faceNorm[i] = dnat[i].reverse().cwiseProduct(Vector2d(1.0, -1.0)) / dnat[i].norm();
    }

    /* Alocate the cell data according to the total cell number */
    data.cell.resize(mesh.totalCells + 1);
    data.tmax = tmax;
    for (std::size_t i = mesh.cellFirstTag; i <= mesh.totalCells; ++i) {
        data.cell[i].RightEV.resize(mesh.dimension);
        data.cell[i].LeftEV.resize(mesh.dimension);
        data.cell[i].face.resize(mesh.numCellFace);
        data.cell[i].qn.resize(3);
    }

    /* Allocate the stencil tag and their vars storage */
    stencils.resize(mesh.stencilWidth + 1);
    qs.resize(NavierStokes.numberOfEquations(), mesh.stencilWidth + 1);
    hs.resize(mesh.stencilWidth + 1);
}

void FiniteVolume::setConstant(const double &SpecificHeat, const double &gasConstant, const double &userCFL) {
    data.gamma = SpecificHeat;
    data.gasConstant = gasConstant;
    data.cfl = userCFL;
    NavierStokes.setConstant(data.gamma, data.gasConstant);
//    initialize(); boundary(); updateVariable(); interpolate(); residual();
//    for (std::size_t i = mesh.cellFirstTag + 1; i <= mesh.cellLastTag; ++i) {
//        std::cout << i << " : " << data.cell[i].residu.transpose() << "\n";
//    }
}

void FiniteVolume::initialize() {
    double pi = 4.0 * atan(1.0);

    if (cases == "Entropy wave") {
        std::cout << "===== Case: Entropy wave =====\n";
        for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
            data.cell[i].pv(rho) = 1.0 + 0.5 * sin(pi * (mesh.cell[i].barycenter(x) +
                                                         mesh.cell[i].barycenter(y)));
            data.cell[i].pv.segment<3>(1) << 1.0, 1.0, 1.0;
            NavierStokes.conservative(data.cell[i].pv, data.cell[i].q);
//            std::cout << data.cell[i].q.transpose() << "\n=====\n";
        }
    } else if (cases == "Bump") {
        std::cout << "===== Case: Bump (Subsonic) =====\n";
        for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
            data.cell[i].pv << 1.4, 1.0, 0.0, 1.0;
            NavierStokes.conservative(data.cell[i].pv, data.cell[i].q);
            NavierStokes.IdealGasEOS(data.cell[i].pv, data.cell[i].T);
            data.cell[i].pvt(0) = 1.5 * pow(data.cell[i].pv.segment<2>(vex).norm() * 0.0, 2);
            data.cell[i].pvt(1) = sqrt(data.cell[i].pvt(0)) * pow(0.09, -0.25) / (0.1 * 3.0);
            kw.conservative(data.cell[i].pvt, data.cell[i].pv(rho), data.cell[i].qt);
        }
    } else if (cases == "Sod") {
        std::cout << "===== Case: Sod Shock Tube =====\n";
        for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
            if (mesh.cell[i].barycenter(x) < 0.5) {
                data.cell[i].pv << 1.0, 0.0, 0.0, 1.0;
            } else {
                data.cell[i].pv << 0.125, 0.0, 0.0, 0.1;
            }
            NavierStokes.conservative(data.cell[i].pv, data.cell[i].q);
        }
    } else if (cases == "Lax") {
        std::cout << "===== Case: Lax =====\n";
        for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
            if (mesh.cell[i].barycenter(x) < 0.5) {
                data.cell[i].pv << 0.445, 0.698, 0.0, 3.528;
            } else {
                data.cell[i].pv << 0.5, 0.0, 0.0, 0.571;
            }
            NavierStokes.conservative(data.cell[i].pv, data.cell[i].q);
        }
    }
    boundary();
}

void FiniteVolume::boundary() {
    std::string boundaryName;
    std::size_t refCell, wallRefCell;
    int refFace;
    double vn;

    for (std::size_t i = mesh.cellLastTag + 1; i <= mesh.totalCells; ++i) {
        boundaryName = mesh.boundaryNames.right.at(mesh.cell[i].bcTag);
        refCell = mesh.cell[i].bcCellRef;
        refFace = mesh.cell[i].bcFaceRef;
        wallRefCell = mesh.cell[i].bcWallRef;

        if (boundaryName == "Inlet") {
            data.cell[i].pv << 1.4, 1.0, 0.0, 1.0;
        } else if (boundaryName == "Outlet" || boundaryName == "ZeroGrad") {
            data.cell[i].pv = data.cell[refCell].pv;
        } else if (boundaryName == "Wall") {
            data.cell[i].pv = data.cell[wallRefCell].pv;
            vn = (data.cell[i].pv(vex) * mesh.cell[refCell].face[refFace].d(y) -
                  data.cell[i].pv(vey) * mesh.cell[refCell].face[refFace].d(x)) /
                 mesh.cell[refCell].face[refFace].ds;
            data.cell[i].pv(vex) = data.cell[i].pv(vex) - 2.0 * vn * mesh.cell[refCell].face[refFace].d(y) /
                                                          mesh.cell[refCell].face[refFace].ds;
            data.cell[i].pv(vey) = data.cell[i].pv(vey) + 2.0 * vn * mesh.cell[refCell].face[refFace].d(x) /
                                                          mesh.cell[refCell].face[refFace].ds;
        }
        NavierStokes.conservative(data.cell[i].pv, data.cell[i].q);
    }
}

void FiniteVolume::updateVariable() {
    Vector2d norm;
    std::size_t refCell;

    for (std::size_t i = mesh.cellFirstTag; i <= mesh.totalCells; ++i) {
        NavierStokes.primitive(data.cell[i].q, data.cell[i].pv);
        NS::Enthalpy(data.cell[i].pv, data.cell[i].q, data.cell[i].enthalpy);
//        std::cout << data.cell[i].enthalpy<< "\n=====\n";
    }

    for (std::size_t i = mesh.cellFirstTag; i <= mesh.totalCells; ++i) {
        for (int f = 1; f < 3; ++f) {
            neighbor = mesh.cell[i].neighbor[f];

            if (neighbor != 0) {
                if (i <= mesh.cellLastTag) {
                    norm = mesh.cell[i].face[f].norm;
                } else {
                    refCell = mesh.cell[i].bcCellRef;
                    norm = mesh.cell[refCell].face[f].norm;
                }

                NavierStokes.RoeAveraging(data.cell[i].enthalpy, data.cell[neighbor].enthalpy, data.cell[i].pv,
                                          data.cell[neighbor].pv, meanVel, meanEnthalpy, meanDynPressure, soundSpd);
                NavierStokes.EigenVectorLeft(meanVel, meanDynPressure, soundSpd, norm,
                                             data.cell[i].LeftEV[f - 1]);
                NavierStokes.EigenVectorRight(meanVel, meanEnthalpy, meanDynPressure, soundSpd,
                                              norm, data.cell[i].RightEV[f - 1]);
//                std::cout << data.cell[i].RightEV[f - 1] << "\n======\n";
            }
        }
    }
}

void FiniteVolume::interpolate() {

    for (std::size_t i = mesh.cellFirstTag; i <= mesh.totalCells; ++i) {
        for (int f = 1; f < 3; ++f) {
            neighbor = mesh.cell[i].neighbor[f];

            if (neighbor != 0) {
                stencils << mesh.cell[i].stencil.segment((f - 1) * mesh.stencilWidth, mesh.stencilWidth),
                        mesh.cell[neighbor].stencil[f * mesh.stencilWidth - 1];
                if (stencils.minCoeff() != 0) {
                    for (int s = 0; s < stencils.size(); ++s) {
                        qs.col(s) = data.cell[i].LeftEV[f - 1] * data.cell[stencils[s]].q;

                        if (stencils[s] <= mesh.cellLastTag) {
                            hs(s) = mesh.cell[stencils[s]].h(f - 1);
                        } else {
                            hs(s) = mesh.cell[mesh.cell[stencils[s]].bcCellRef].h(f - 1);
                        }
                    }

                    for (int eq = 0; eq < NavierStokes.numberOfEquations(); ++eq) {
                        schemes.WENOJS_5th(qs.row(eq),
                                           data.cell[i].face[f].ql(eq), data.cell[i].face[f].qr(eq));
                    }

                    data.cell[i].face[f].ql = data.cell[i].RightEV[f - 1] * data.cell[i].face[f].ql;
                    data.cell[i].face[f].qr = data.cell[i].RightEV[f - 1] * data.cell[i].face[f].qr;
                    NavierStokes.primitive(data.cell[i].face[f].ql, pvl);
                    NavierStokes.primitive(data.cell[i].face[f].qr, pvr);

                    for (int dir = x; dir <= y; ++dir) {
                        data.cell[i].face[f].flux.col(dir) = NavierStokes.RusanovFlux(data.cell[i].face[f].ql, pvl,
                                                                                      data.cell[i].face[f].qr, pvr,
                                                                                      data.cell[i].LeftEV[f - 1],
                                                                                      data.cell[i].RightEV[f - 1],
                                                                                      dir);
                    }
                }
            }
        }
    }
}

void FiniteVolume::residual() {
    boundary();
    updateVariable();
    interpolate();

    for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
        for (int f = 1; f < 3; ++f) {
            neighbor = mesh.cell[i].neighbor[circularNumber[f + 2]];
            data.cell[i].face[circularNumber[f + 2]].flux = data.cell[neighbor].face[f].flux;
        }
    }

    for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
        for (int f = 0; f < mesh.numCellFace; ++f) {
            res.col(f) = data.cell[i].face[f].flux.col(x) * mesh.cell[i].face[f].d(y) -
                         data.cell[i].face[f].flux.col(y) * mesh.cell[i].face[f].d(x);
        }
        data.cell[i].residu = -res.rowwise().sum() / mesh.cell[i].area;
//        std::cout << data.cell[i].residu.transpose() << "\n======\n";
    }
}

void FiniteVolume::timeStepping() {
    double scale = 3.0 / 3.0;

    while (data.t < data.tmax) {
        globalMaxWaveSpd.setZero();
        for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
            waveSpd = sqrt(fabs(data.gamma * data.cell[i].pv(p) / data.cell[i].pv(rho)));

            for (int dir = x; dir <= y; ++dir) {
                globalMaxWaveSpd(dir) = fmax(globalMaxWaveSpd(dir), fmax(fabs(data.cell[i].pv(dir + 1) - waveSpd),
                                                                         fabs(data.cell[i].pv(dir + 1) + waveSpd)));
            }
        }

        data.dt = data.cfl * pow(mesh.minGridSize, scale) / globalMaxWaveSpd.maxCoeff();
        data.dt = data.t + data.dt > data.tmax ? data.tmax - data.t : data.dt;

        stage = 1;
        residual();
        for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
            data.cell[i].qn[0] = data.cell[i].q;
            data.cell[i].qn[1] = ssprk.RK_3_3<Vector4d>(data.cell[i].qn[0], data.cell[i].qn[0],
                                                        data.cell[i].residu, data.dt, stage);
            data.cell[i].q = data.cell[i].qn[1];
        }

        stage += 1;
        residual();
        for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
            data.cell[i].qn[2] = ssprk.RK_3_3<Vector4d>(data.cell[i].qn[0], data.cell[i].qn[1],
                                                        data.cell[i].residu, data.dt, stage);
            data.cell[i].q = data.cell[i].qn[2];
        }

        stage += 1;
        residual();
        for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
            data.cell[i].q = ssprk.RK_3_3<Vector4d>(data.cell[i].qn[0], data.cell[i].qn[2],
                                                    data.cell[i].residu, data.dt, stage);
        }

        data.t += data.dt;
        cout << data.t << "\n";
    }

    for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
        NavierStokes.primitive(data.cell[i].q, data.cell[i].pv);
    }
    double err = 0.0;
    double pi = 4.0 * atan(1.0);
    for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
        err += fabs(data.cell[i].pv(rho) - (1.0 + 0.5 * sin(pi * (mesh.cell[i].barycenter(x) - data.tmax +
                                                                  mesh.cell[i].barycenter(y) - data.tmax)))) /
               mesh.internalCells;
    }
    std::cout << "err: " << err << "\n";
}

void FiniteVolume::run() {
    initialize();
    timeStepping();
}

void FiniteVolume::save(const std::string &filename) {
    std::ofstream ofile;
    const char separator = ' ';
    const int strwidth = 20;
    const int numwidth = 20;
    const int pn = 8;

    ofile.open(filename);

    if (ofile.good()) {
        ofile << "t: " << data.t << "\n";
        ofile << std::right << std::setw(strwidth) << std::setfill(separator) << "x";
        ofile << std::right << std::setw(strwidth) << std::setfill(separator) << "y";
        ofile << std::right << std::setw(strwidth) << std::setfill(separator) << "Density";
        ofile << std::right << std::setw(strwidth) << std::setfill(separator) << "u";
        ofile << std::right << std::setw(strwidth) << std::setfill(separator) << "v";
        ofile << std::right << std::setw(strwidth) << std::setfill(separator) << "Pressure";
        ofile << "\n";
    } else {
        exit(1);
    }

    for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
        if (round(mesh.cell[i].barycenter(y) * 10000.0) == 1525.0) {
            ofile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
                  << mesh.cell[i].barycenter[x];
            ofile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
                  << mesh.cell[i].barycenter[y];
            ofile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
                  << data.cell[i].pv(rho);
            ofile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
                  << data.cell[i].pv(vex);
            ofile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
                  << data.cell[i].pv(vey);
            ofile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
                  << data.cell[i].pv(p);
            ofile << "\n";
        }
    }
}

FiniteVolume::~FiniteVolume() {
    std::cout << "Destructor of finite volume class has been called.\n";
}

//
// Created by Indra on 02/10/20.
//

#ifndef FREEFLOW_K_OMEGA_H
#define FREEFLOW_K_OMEGA_H

#include "../include/EigenType.h"

class k_Omega {
public:
    k_Omega();
    void conservative(const Vector2d &pvt, const double &density, Vector2d &qt);
    friend double KroneckerDelta(const int &i, const int &j);

private:
    double Pr, Prt, C1, C2;
};


#endif //FREEFLOW_K_OMEGA_H

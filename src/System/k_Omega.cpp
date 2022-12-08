//
// Created by Indra on 02/10/20.
//

#include "k_Omega.h"

k_Omega::k_Omega() {
    Pr = 0.72, Prt = 0.9;
    C1 = 1.458E-6, C2 = 110.4;
}

double KroneckerDelta(const int &i, const int &j) {
    return i != j ? 0.0 : 1.0;
}

void k_Omega::conservative(const Vector2d &pvt, const double &density, Vector2d &qt) {
    qt = pvt * density;
}

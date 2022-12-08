//
// Created by mf on 29/09/20.
//

#include "SSPRK.h"

template<typename T>
T SSPRK::RK_3_3(const T &q0, const T &qn, const T &rhs, const double &dt, const int &stage) {
    if (stage == 1) {
        return q0 + dt * rhs;
    } else if (stage == 2) {
        return 0.75 * q0 + 0.25 * (qn + rhs * dt);
    } else {
        return (q0 + 2.0 * (qn + dt * rhs)) / 3.0;
    }
}

template Vector4d SSPRK::RK_3_3<Vector4d>(const Vector4d &, const Vector4d &, const Vector4d &,
                                          const double &, const int &);
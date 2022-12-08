//
// Created by mf on 29/09/20.
//

#ifndef FREEFLOW_SSPRK_H
#define FREEFLOW_SSPRK_H

#include "../include/EigenType.h"
#include <vector>

class SSPRK {
public:
    template <typename T> T RK_3_3(const T& q0, const T& qn, const T& rhs,
                                   const double& dt, const int& stage);
};

#endif //FREEFLOW_SSPRK_H

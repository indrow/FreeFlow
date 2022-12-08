//
// Created by Indra on 8/6/2020.
//

#ifndef FREEFLOW_EIGENTYPE_H
#define FREEFLOW_EIGENTYPE_H

#include <Eigen/Dense>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
typedef Eigen::Matrix<double, 5, 1> Vector5d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 9, 1> Vector9d;
typedef Eigen::Matrix4d Matrix4d;

class hasher {
public:
    std::size_t operator()(std::vector<std::size_t> const &vec) const {
        std::size_t seed = vec.size();
        for (auto &i : vec) {
            seed ^= i + 0x9e37779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};
#endif //FREEFLOW_EIGENTYPE_H

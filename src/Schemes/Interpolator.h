//
// Created by Indra on 01/10/20.
//

#ifndef FREEFLOW_INTERPOLATOR_H
#define FREEFLOW_INTERPOLATOR_H

#include "../include/EigenType.h"
#include <iostream>

class Interpolator {
public:
    Interpolator();

    void WENOJS_3rd(const Vector4d &q, double &ql, double &qr);

    void WENONP_3rd(const Vector4d &q, double &ql, double &qr);

    void WENOZQ_3rd(const Vector4d &q, double &ql, double &qr);

    void WENOJS_5th(const Vector6d &q, double &ql, double &qr);

    void WENOZQ_5th(const Vector6d &q, double &ql, double &qr);

    void WENONM_5th(const Vector6d &q, const Vector6d &h, double &ql, double &qr);

    ~Interpolator();

private:
    const double eps, epsz;

    double tau;
    Matrix2d cplus_3rd, cmin_3rd;
    Vector2d dplus_3rd, dmin_3rd, p_3rd, si_3rd, weights_3rd;
    Vector3d dq_3rd, pq_3rd, wq_3rd;

    Vector2d SI_3rd(const Vector3d &q);
    double SI_fs_3rd(const Vector3d &q) const;
    static double cutoff(const double &val, const double &ct);

    double csi_5th;
    Matrix3d cplus_5th, cmin_5th;
    Vector3d dplus_5th, dmin_5th, p_5th, si_5th, weights_5th;
    Vector5d dq_5th, pq_5th, wq_5th, cqmin_5th, cqplus_5th;

    Vector3d SI_5th(const Vector5d &q);
    static double SI_fs_5rd(const Vector5d &q) ;

    double alpha;
    Vector3d qm, qmt;
    Matrix3d cplus_nm_5th, cmin_nm_5th;

    double linear_interpolation(const double &q1, const double &q2,
                                const double &h1, const double &h2);
};


#endif //FREEFLOW_INTERPOLATOR_H

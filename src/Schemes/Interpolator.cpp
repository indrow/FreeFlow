//
// Created by Indra on 01/10/20.
//

#include "Interpolator.h"

/*
 * Sign convensions :
 * qr = q_{i+1/2}^{+}
 * ql  = q_{i+1/2)^{-}
 * */

Interpolator::Interpolator() : eps(1.0E-6), epsz(1.0E-40) {
    std::cout << "Constructor of the Interpolator has been called.\n";
    tau = 0.0;

    cmin_3rd << -0.5, 1.5, 0.5, 0.5;
    cplus_3rd << 0.5, 0.5, 1.5, -0.5;
    dmin_3rd << 1.0 / 3.0, 2.0 / 3.0;
    dplus_3rd << 2.0 / 3.0, 1.0 / 3.0;
    dq_3rd << 0.92, 0.04, 0.04;

    cmin_5th << 1.0 / 3.0, -7.0 / 6.0, 11.0 / 6.0,
            -1.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0,
            1.0 / 3.0, 5.0 / 6.0, -1.0 / 6.0;
    cplus_5th = cmin_5th;
    cplus_5th.row(0).swap(cplus_5th.row(2));
    cplus_5th.col(0).swap(cplus_5th.col(2));
    cqmin_5th << 1.0 / 30.0, -13.0 / 60.0, 47.0 / 60.0, 9.0 / 20.0, -1.0 / 20.0;
    cqplus_5th << -1.0 / 20.0, 9.0 / 20.0, 47.0 / 60.0, -13.0 / 60.0, 1.0 / 30.0;
    dq_5th << 0.86, 0.02, 0.02, 0.05, 0.05;

    dmin_5th << 0.1, 0.6, 0.3;
    dplus_5th << 0.3, 0.6, 0.1;

    csi_5th = 13.0 / 12.0;

    alpha = 0.5;
    cmin_nm_5th << 2.0, -1.0 / 3.0, -2.0 / 3.0,
            -1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0,
            2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0;
    cplus_nm_5th = cmin_nm_5th;
    cplus_nm_5th.row(0).swap(cplus_nm_5th.row(2));
    cplus_nm_5th.col(0).segment(0, 2).swap(cplus_nm_5th.col(2).segment(0, 2));
}

void Interpolator::WENOJS_3rd(const Vector4d &q, double &ql, double &qr) {
    p_3rd(0) = cmin_3rd.row(0).dot(q.segment(0, 2));
    p_3rd(1) = cmin_3rd.row(1).dot(q.segment(1, 2));

    weights_3rd = dmin_3rd.array() / (eps + SI_3rd(q.segment(0, 3)).array()).pow(2);
    weights_3rd = weights_3rd.array() / weights_3rd.sum();
    ql = p_3rd.dot(weights_3rd);
    ql = std::isnan(ql) ? 0.0 : ql;

    p_3rd(0) = cplus_3rd.row(0).dot(q.segment(1, 2));
    p_3rd(1) = cplus_3rd.row(1).dot(q.segment(2, 2));

    weights_3rd = dplus_3rd.array() / (eps + SI_3rd(q.segment(1, 3)).array()).pow(2);
    weights_3rd = weights_3rd.array() / weights_3rd.sum();
    qr = p_3rd.dot(weights_3rd);
    qr = std::isnan(ql) ? 0.0 : qr;
}

void Interpolator::WENONP_3rd(const Vector4d &q, double &ql, double &qr) {
    Vector2d si_3rdl;

    p_3rd(0) = cmin_3rd.row(0).dot(q.segment(0, 2));
    p_3rd(1) = cmin_3rd.row(1).dot(q.segment(1, 2));

    si_3rdl = SI_3rd(q.segment(0, 3));
    tau = fabs(si_3rdl.sum() / 2.0 - SI_fs_3rd(q.segment(0, 3)));

    weights_3rd = (1.0 + tau / (si_3rdl.array() + epsz).pow(0.75));
    weights_3rd = dmin_3rd.cwiseProduct(weights_3rd);
    weights_3rd = weights_3rd.array() / weights_3rd.sum();
    ql = p_3rd.dot(weights_3rd);
    ql = std::isnan(ql) ? 0.0 : ql;

    p_3rd(0) = cplus_3rd.row(0).dot(q.segment(1, 2));
    p_3rd(1) = cplus_3rd.row(1).dot(q.segment(2, 2));

    si_3rdl = SI_3rd(q.segment(1, 3));
    tau = fabs(si_3rdl.sum() / 2.0 - SI_fs_3rd(q.segment(1, 3)));

    weights_3rd = (1.0 + tau / (si_3rdl.array() + epsz).pow(0.75));
    weights_3rd = dplus_3rd.cwiseProduct(weights_3rd);
    weights_3rd = weights_3rd.array() / weights_3rd.sum();
    qr = p_3rd.dot(weights_3rd);
    qr = std::isnan(qr) ? 0.0 : qr;
}

void Interpolator::WENOZQ_3rd(const Vector4d &q, double &ql, double &qr) {
    Vector3d si_3rdl;

    pq_3rd(0) = cmin_5th.row(1).dot(q.segment(0, 3));
    pq_3rd(1) = cmin_3rd.row(0).dot(q.segment(0, 2));
    pq_3rd(2) = cmin_3rd.row(1).dot(q.segment(1, 2));

    si_3rdl(0) = SI_fs_3rd(q.segment(0, 3));
    si_3rdl.segment<2>(1) = SI_3rd(q.segment(0, 3));
//    tau = fabs(si_3rdl(0) - si_3rdl.segment<2>(1).sum() / 2.0);
    tau = (fabs(si_3rdl.segment<2>(1).sum() / 2.0 - 0.25 * pow(q(2) - q(0), 2)));

    wq_3rd = dq_3rd.array() * (1.0 + tau / (si_3rdl.array() + epsz));
    wq_3rd = wq_3rd.array() / wq_3rd.sum();
    ql = wq_3rd(0) / dq_3rd(0) * (pq_3rd(0) - dq_3rd.segment<2>(1).dot(pq_3rd.segment<2>(1))) +
         wq_3rd.segment<2>(1).dot(pq_3rd.segment<2>(1));

    pq_3rd(0) = cplus_5th.row(1).dot(q.segment(1, 3));
    pq_3rd(1) = cplus_3rd.row(0).dot(q.segment(1, 2));
    pq_3rd(2) = cplus_3rd.row(1).dot(q.segment(2, 2));

    si_3rdl(0) = SI_fs_3rd(q.segment(1, 3));
    si_3rdl.segment<2>(1) = SI_3rd(q.segment(1, 3));
//    tau = fabs(si_3rdl(0) - si_3rdl.segment<2>(1).sum() / 2.0);
    tau = (fabs(si_3rdl.segment<2>(1).sum() / 2.0 - 0.25 * pow(q(3) - q(1), 2)));

    wq_3rd = dq_3rd.array() * (1.0 + tau / (si_3rdl.array() + epsz));
    wq_3rd = wq_3rd.array() / wq_3rd.sum();
    qr = wq_3rd(0) / dq_3rd(0) * (pq_3rd(0) - dq_3rd.segment<2>(1).dot(pq_3rd.segment<2>(1))) +
         wq_3rd.segment<2>(1).dot(pq_3rd.segment<2>(1));
}

Vector2d Interpolator::SI_3rd(const Vector3d &q) {
    si_3rd(0) = q(0) - q(1);
    si_3rd(1) = q(1) - q(2);

    return si_3rd.array().pow(2);
}

double Interpolator::SI_fs_3rd(const Vector3d &q) const {
    return csi_5th * pow(q(0) - 2.0 * q(1) + q(2), 2) + 0.25 * pow(q(0) - q(2), 2);
}

double Interpolator::cutoff(const double &val, const double &ct) {
    return val < ct ? 0.0 : 1.0;
}

void Interpolator::WENOJS_5th(const Vector6d &q, double &ql, double &qr) {
    p_5th(0) = cmin_5th.row(0).dot(q.segment(0, 3));
    p_5th(1) = cmin_5th.row(1).dot(q.segment(1, 3));
    p_5th(2) = cmin_5th.row(2).dot(q.segment(2, 3));

    weights_5th = dmin_5th.array() / (eps + SI_5th(q.segment(0, 5)).array()).pow(2);
    weights_5th = weights_5th.array() / weights_5th.sum();
    ql = p_5th.dot(weights_5th);

    p_5th(0) = cplus_5th.row(0).dot(q.segment(1, 3));
    p_5th(1) = cplus_5th.row(1).dot(q.segment(2, 3));
    p_5th(2) = cplus_5th.row(2).dot(q.segment(3, 3));

    weights_5th = dplus_5th.array() / (eps + SI_5th(q.segment(1, 5)).array()).pow(2);
    weights_5th = weights_5th.array() / weights_5th.sum();
    qr = p_5th.dot(weights_5th);
}



void Interpolator::WENONM_5th(const Vector6d &q, const Vector6d &h, double &ql, double &qr) {
    qm(0) = q(2);
    qm(1) = linear_interpolation(q(1), q(2), h(1), h(2));
    qm(2) = ((h(0) + 2.0 * h(1)) * q(1) - h(1) * q(0)) / (h(0) + h(1));
    p_5th(0) = cmin_nm_5th.row(0).dot(qm);
    si_5th(0) = pow(2.0 * qm(0) - qm(1) - qm(2), 2) + csi_5th * pow(2.0 * (qm(1) - qm(2)), 2);
    qmt(0) = qm(1);

    qm(0) = qm(1);
    qm(1) = q(2);
    qm(2) = linear_interpolation(q(2), q(3), h(2), h(3));
    p_5th(1) = cmin_nm_5th.row(1).dot(qm);
    si_5th(1) = pow(-qm(0) + qm(2), 2) + csi_5th * pow(2.0 * (qm(0) + qm(2)) - 4.0 * qm(1), 2);
    qmt(1) = qm(2);

    qm(0) = qm(2);
    qm(1) = q(3);
    qm(2) = linear_interpolation(q(3), q(4), h(3), h(4));
    p_5th(2) = cmin_nm_5th.row(2).dot(qm);
    si_5th(2) = pow(-3.0 * qm(0) + 4.0 * qm(1) - qm(2), 2) + csi_5th * pow(2.0 * (qm(0) + qm(2)) - 4.0 * qm(1), 2);
    qmt(2) = qm(2);

    weights_5th = dmin_5th.array() / (eps + si_5th.array()).pow(2);
    weights_5th = weights_5th.array() / weights_5th.sum();
    ql = p_5th.dot(weights_5th);

    qm(0) = qmt(0);
    qm(1) = q(2);
    qm(2) = qmt(1);
    p_5th(0) = cplus_nm_5th.row(0).dot(qm);
    si_5th(0) = pow(qm(0) - 4.0 * qm(1) + 3.0 * qm(2), 2) + csi_5th * pow(2.0 * (qm(0) + qm(2)) - 4.0 * qm(1), 2);

    qm(0) = qmt(1);
    qm(1) = q(3);
    qm(2) = qmt(2);
    p_5th(1) = cplus_nm_5th.row(1).dot(qm);
    si_5th(1) = pow(-qm(0) + qm(2), 2) + csi_5th * pow(2.0 * (qm(0) + qm(2)) - 4.0 * qm(1), 2);

    qm(0) = q(3);
    qm(1) = qmt(2);
    qm(2) = ((h(5) + 2.0 * h(4)) * q(4) - h(4) * q(5)) / (h(4) + h(5));
    p_5th(2) = cplus_nm_5th.row(2).dot(qm);
    si_5th(2) = pow(2.0 * qm(0) - qm(1) - qm(2), 2) + csi_5th * pow(2.0 * (qm(1) - qm(2)), 2);
    weights_5th = dplus_5th.array() / (eps + si_5th.array()).pow(2);
    weights_5th = weights_5th.array() / weights_5th.sum();
    qr = p_5th.dot(weights_5th);
}

void Interpolator::WENOZQ_5th(const Vector6d &q, double &ql, double &qr) {
    Vector3d si_thjs;
    Vector5d si_5thl, taul;

    pq_5th(0) = cqmin_5th.dot(q.segment<5>(0));
    pq_5th(1) = cmin_3rd.row(0).dot(q.segment<2>(1));
    pq_5th(2) = cmin_3rd.row(1).dot(q.segment<2>(2));
    pq_5th(3) = cmin_5th.row(0).dot(q.segment<3>(0));
    pq_5th(4) = cmin_5th.row(2).dot(q.segment<3>(2));
    si_thjs = SI_5th(q.segment<5>(0));

    si_5thl(0) = SI_fs_5rd(q.segment<5>(0));
    si_5thl.segment<2>(1) = SI_3rd(q.segment<3>(1));
    si_5thl.segment<2>(3) << si_thjs(0), si_thjs(2);

    tau = fabs(si_5thl(0) - (si_5thl(3) + si_5thl(4) + 4.0 * si_thjs(1)) / 6.0);
    taul(0) = tau, taul(3) = tau, taul(4) = tau;
    tau = pow((fabs(si_5thl(0) - si_5thl(1)) + fabs(si_5thl(0) - si_5thl(2))) / 2.0, 2);
    taul(1) = tau, taul(2) = tau;

    wq_5th = dq_5th.array() * (1.0 + taul.array() / (si_5thl.array() + epsz));
    wq_5th = wq_5th.array() / wq_5th.sum();
//    wq_5th = wq_5th.array() / wq_5th.segment<3>(0).sum();
    ql = wq_5th(0) / dq_5th(0) * (pq_5th(0) - pq_5th.segment<4>(1).dot(dq_5th.segment<4>(1))) +
         pq_5th.segment<4>(1).dot(wq_5th.segment<4>(1));

    pq_5th(0) = cqplus_5th.dot(q.segment<5>(1));
    pq_5th(1) = cplus_3rd.row(0).dot(q.segment<2>(2));
    pq_5th(2) = cplus_3rd.row(1).dot(q.segment<2>(3));
    pq_5th(3) = cplus_5th.row(0).dot(q.segment<3>(1));
    pq_5th(4) = cplus_5th.row(2).dot(q.segment<3>(3));
    si_thjs = SI_5th(q.segment<5>(1));

    si_5thl(0) = SI_fs_5rd(q.segment<5>(1));
    si_5thl.segment<2>(1) = SI_3rd(q.segment<3>(2));
    si_5thl.segment<2>(3) << si_thjs(0), si_thjs(2);

    tau = fabs(si_5thl(0) - (si_5thl(3) + si_5thl(4) + 4.0 * si_thjs(1)) / 6.0);
    taul(0) = tau, taul(3) = tau, taul(4) = tau;
    tau = pow((fabs(si_5thl(0) - si_5thl(1)) + fabs(si_5thl(0) - si_5thl(2))) / 2.0, 2);
    taul(1) = tau, taul(2) = tau;

    wq_5th = dq_5th.array() * (1.0 + taul.array() / (si_5thl.array() + epsz));
    wq_5th = wq_5th.array() / wq_5th.sum();
//    wq_5th = wq_5th.array() / wq_5th.segment<3>(0).sum();

    qr = wq_5th(0) / dq_5th(0) * (pq_5th(0) - pq_5th.segment<4>(1).dot(dq_5th.segment<4>(1))) +
         pq_5th.segment<4>(1).dot(wq_5th.segment<4>(1));
}

Vector3d Interpolator::SI_5th(const Vector5d &q) {
    si_5th(0) = pow(0.5 * (q(0) - 4.0 * q(1) + 3.0 * q(2)), 2) +
                csi_5th * pow(q(0) - 2.0 * q(1) + q(2), 2);
    si_5th(1) = pow(0.5 * (-q(1) + q(3)), 2) + csi_5th * pow(q(1) - 2.0 * q(2) + q(3), 2);
    si_5th(2) = pow(0.5 * (-3.0 * q(2) + 4.0 * q(3) - q(4)), 2) +
                csi_5th * pow(q(2) - 2.0 * q(3) + q(4), 2);

    return si_5th;
}

double Interpolator::SI_fs_5rd(const Vector5d &q) {
    return (q(4) * (6908.0 * q(4) - 51001.0 * q(3) + 67923.0 * q(2) - 38947.0 * q(1) + 8209.0 * q(0)) +
            q(3) * (104963.0 * q(3) - 299076.0 * q(2) + 179098.0 * q(1) - 38947.0 * q(0)) +
            q(2) * (231153.0 * q(2) - 299076.0 * q(1) + 67923.0 * q(0)) +
            q(1) * (104963.0 * q(1) - 51001.0 * q(0)) + 6908.0 * q(0) * q(0)) / 5040.0;
}

double Interpolator::linear_interpolation(const double &q1, const double &q2,
                                          const double &h1, const double &h2) {
    alpha = h1 / (h1 + h2);
    return (1.0 - alpha) * q1 + alpha * q2;
}

Interpolator::~Interpolator() {
    std::cout << "Destructor of Interpolator has been called.\n";
}
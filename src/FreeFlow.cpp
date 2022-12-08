#include "Solver/FiniteVolume.h"

int main() {
    std::string meshfile = "../UserInput/test_periodic_50.msh",
                cases = "Entropy wave", file = "Lax_WENOZQ3_1D.dat";
    double gamma = 1.4, cfl = 0.4, tmax = 1.0;
    double gasConstant = 8.3145;
    int stn = 5, numPeriodicBc = 2;

    FiniteVolume fv(meshfile, cases, stn, numPeriodicBc, tmax);
    fv.setConstant(gamma, gasConstant, cfl);
//    fv.initialize();
//    fv.boundary();
//    fv.updateVariable();
//    fv.interpolate();
//    fv.residual();
    fv.run();
//    fv.save(file);

    /* Testing */
//    Interpolator ip;
//    double ql, qr;
//    Vector6d q, h;
//    q << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
//    h << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
//    ip.WENOZQ_5th(q, ql, qr);
//    std::cout << ql << "  " << qr << std::endl;

    return 0;
}

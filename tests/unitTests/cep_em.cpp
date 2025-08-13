#include "celec_mech.h"
#include "ComMod.h"
#include "test_common.h"
#include <cmath>
#include <math.h>

class EMTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(EMTest, LandModel) {
    // Test implementation for LandModel
    celec_mech mech;

    double dlambda_dt = 0.1; 
    double nominal_lambda = 1; 
    double T = 0.0;
    double Ta = 0.0;
    double Tp = 0.0;

    int const nY = 7; 
    Vector<double> Y0(nY);
    Y0(0) = 0.0052; 
    Y0(1) = 0.0444; 
    Y0(2) = 0.2609; 
    Y0(3) = 0.8900; 
    Y0(4) = 0.0000;
    Y0(5) = 0.0000;
    Y0(6) = 0.0000; // Assuming Y(6) is also needed

    std::vector<double> dYactual = {4.3920e-04, 0.0029, 0.0029, -0.0066, 1.0, 1.0, 0.0}; 

    Vector<double> dY(nY);
    dY(0) = 0.0;
    dY(1) = 0.0;
    dY(2) = 0.0;
    dY(3) = 0.0;
    dY(4) = 0.0;
    dY(5) = 0.0;
    dY(6) = 0.0;

    double Cai = 0.5040; 

    mech.land_model(Y0, dY, Cai, nominal_lambda, dlambda_dt, T, Ta, Tp);

    // cout each component of Y0
    for (int i = 0; i < nY; ++i) {
        std::cout << "dY(" << i << ") = " << dY(i) << " dYactual(" << i << ") = " << dYactual[i] << std::endl;
    }

}

TEST_F(EMTest, SingleMyocyteContraction) {
    // initialize state
    int const nY = 7; 
    Vector<double> Y0(nY);
    celec_mech mech;

    Y0(0) = 0.0;  // Initial value for Y(0)
    Y0(1) = 0.0;  // Initial value for Y(1)
    Y0(2) = 0.0;  // Initial value for Y(2)
    Y0(3) = 0.0;  // Initial value for Y(3)
    Y0(4) = 0.0;  // Initial value for Y(4)
    Y0(5) = 0.0;  // Initial value for Y(5)
    Y0(6) = 0.0;  // Initial value for Y(6)

    Vector<double> Y = Y0;
    double F = 0.0, Ta = 0.0, Tp = 0.0;
    double T; 
    int const nsteps = 1000; 
    double t = 0.0;
    double dt = (M_PI/5)/nsteps; 

    std::vector<double> Factual = {1.2594e-04, 6.5044e-04, 0.0018, 0.0047, 0.0126,
        0.0271, 0.0467, 0.0675, 0.0855, 0.0980, 0.1039, 0.1028, 0.0947,
        0.0802, 0.0611, 0.0407, 0.0232, 0.0119, 0.0059, 0.0027};

    for (size_t step = 0; step < nsteps; ++step) {
        double Cai = (1/(2.5)) * sin(5*t) + 1/6;
        // if dlambda_dt or lambda vary, you can compute them here; use constants for now
        double lambda = 1.0;
        double dlambda_dt = 0.0;
        mech.integ_rk(nY, Y, T, Ta, Tp, dt, Cai, lambda, dlambda_dt);
        t += dt;
        if (step % 50 == 0) {
            // print out
            std::cout << "Step " << step << ": Y(0) = " << Y(0) << ", Factual = " << Factual[step / 50] << std::endl;
            // std::cout << "Step " << step << ": Y(1) = " << Y(1) << std::endl;
            // std::cout << "t = " << t << std::endl;
            // ASSERT_NEAR(Y(0), Factual[step / 50], 1e-2) << "Failed at step " << step;
        }
    }
}
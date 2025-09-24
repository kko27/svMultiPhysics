/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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
    int const nsteps = 1001; 
    double t = 0.0; 
    double tCai = 0.0; 
    double dt = 1.0; 
    double dtCai = (M_PI/5); 

    std::vector<double> Factual = {7.0601e-04, 0.004155, 0.002316, 0.002557, 0.005096, 0.012981, 0.027499,
        0.047117, 0.067922, 0.085835, 0.098221, 0.103934, 0.102656,
        0.094444, 0.079856, 0.060638, 0.040262, 0.022969, 0.011724,
        0.005796, 0.002701};

    for (size_t step = 0; step < nsteps; ++step) {
        double Cai = (1.0/2.5) * sin(5*tCai/1000.0) + 1.0/6;
        // if dlambda_dt or lambda vary, you can compute them here; use constants for now
        double lambda = 1.0;
        double dlambda_dt = 0.0;
        mech.integ_rk(nY, Y, T, Ta, Tp, dt, Cai, lambda, dlambda_dt);
        t += dt;
        tCai += dtCai;
        if (step % 50 == 0) {
            // print out
            // std::cout << "time " << t << std::endl;
            // std::cout << "tCai: " << tCai << std::endl;
            // std::cout << "Step " << step << ": Y(0) = " << Y(0) << std::endl;
            // std::cout << "Factual = " << Factual[step / 50] << std::endl;
            ASSERT_NEAR(Y(0), Factual[step / 50], 1e-3) << "Failed at step " << step;
        }
    }
}
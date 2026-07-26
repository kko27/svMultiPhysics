/* Copyright (c) Stanford University, The Regents of the University of
 * California, and others.
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

#include "../test_common.h"
#include "active_stress_land_niederer.h"

class TestableLandNiederer : public LandNiederer {
public:
  Vector<double> evaluate_getf(const double t, const Vector<double> &state,
                               const double calcium,
                               const double fiber_stretch,
                               const double fiber_stretch_rate) const {
    return getf(t, state, calcium, fiber_stretch, fiber_stretch_rate);
  }

  void set_parameters_for_test() {
    CaRef = 0.805;
    eta_Tm = 5.0;
    k_uw = 0.182;
    k_ws = 0.012;
    Tref = 120.0;
    k_TRPN = 0.1;
    eta_TRPN = 2.0;
    k_u = 1.0;
    TRPN50 = 0.35;
    rw = 0.5;
    rs = 0.25;
    gamma_s = 0.0085;
    gamma_w = 0.615;
    phi = 2.23;
    Aeff = 25.0;
    beta0 = 2.3;
    beta1 = -2.4;
  }
};

class ElectromechanicsModelTest : public ::testing::Test {
protected:
  void SetUp() override { model.set_parameters_for_test(); }

  void TearDown() override {}

  static Vector<double> initial_state() {
    return {0.8900, 0.0444, 0.2609, 0.0052, 0.0000, 0.0000};
  }

  static Vector<double> expected_rhs() {
    return {-0.0066, 0.0029, 0.0029, 0.0004, 1.0, 1.0};
  }

  TestableLandNiederer model;
};

TEST_F(ElectromechanicsModelTest, LandNiedererGetf) {
  // Evaluate the Land-Niederer ODE right-hand side at one state and compare
  // against reference values computed from the MATLAB implementation.

  const double t = 0.0;
  const double calcium = 0.5040;
  const double fiber_stretch = 1.0;
  const double fiber_stretch_rate = 0.1;

  const Vector<double> state = initial_state();
  const Vector<double> expected = expected_rhs();

  const Vector<double> rhs =
      model.evaluate_getf(t, state, calcium, fiber_stretch,
                          fiber_stretch_rate);

  ASSERT_EQ(rhs.size(), expected.size());

  for (int i = 0; i < rhs.size(); ++i) {
    ASSERT_NEAR(rhs[i], expected[i], 1e-3)
        << "RHS mismatch at component " << i;
  }
}

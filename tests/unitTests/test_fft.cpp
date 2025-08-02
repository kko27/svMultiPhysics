#include "fft.h"
#include "ComMod.h"
#include <gtest/gtest.h>
#include <cmath>

class FFTTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(FFTTest, LinearRamp) {
    // Example: 1D signal, linear ramp: y = 2*t + 1
    // t: 0, 1, 2, 3
    // y: 1, 3, 5, 7
    std::vector<std::vector<double>> temporal_values = {
        {0.0, 1.0},
        {1.0, 3.0},
        {2.0, 5.0},
        {3.0, 7.0}
    };
    int np = temporal_values.size();

    fcType gt;
    gt.d = 1; // 1D signal
    gt.n = 2; // Compute 2 Fourier coefficients
    gt.qi.resize(gt.d);
    gt.qs.resize(gt.d);
    gt.r.resize(gt.d, gt.n);
    gt.i.resize(gt.d, gt.n);

    fft(np, temporal_values, gt);

    // For a linear ramp y = 2*t + 1, the DC component (n=0) should be the mean value
    // Over t in [0,3]: mean = (1+3+5+7)/4 = 4
    // The slope is 2, so qs should be 2
    ASSERT_NEAR(gt.qi[0], 1.0, 1e-10) << "gt.qi[0] should be 1.0 (intercept)";
    ASSERT_NEAR(gt.qs[0], 2.0, 1e-10) << "gt.qs[0] should be 2.0 (slope)";
    // The DC Fourier coefficient (real part, n=0) should be close to 0 (since the ramp is removed)
    ASSERT_NEAR(gt.r(0,0), 0.0, 1e-10) << "gt.r(0,0) should be ~0 after ramp removal";
    // The imaginary part should be zero for real signals
    ASSERT_NEAR(gt.i(0,0), 0.0, 1e-10) << "gt.i(0,0) should be 0";
}

// main() is provided by gtest 
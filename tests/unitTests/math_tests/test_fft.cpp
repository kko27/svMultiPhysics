#include "fft.h"
#include "ComMod.h"
#include "../test_common.h"
#include <cmath>

class FFTTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(FFTTest, SinCosLinearCombination) {
    int N = 100;
    double x_start = 0.0;
    double x_end = 10.0;

    std::vector<std::vector<double>> temporal_values;
    temporal_values.reserve(N);

    double step = (x_end - x_start) / (N - 1);
    for (int i = 0; i < N; ++i) {
        double t = x_start + i * step;
        double y = std::sin(t) + std::cos(t) + 0.1 * t;
        temporal_values.push_back({t, y});
    }
    fcType gt;
    gt.d = 1;
    gt.n = 16;
    gt.qi.resize(gt.d);
    gt.qs.resize(gt.d);
    gt.r.resize(gt.d, gt.n);
    gt.i.resize(gt.d, gt.n);

    fft(N, temporal_values, gt);

    ASSERT_NEAR(gt.qs[0], -0.13830, 1e-2) << "Expected slope ~-0.13830";
    // Check the first three Fourier coefficients
    ASSERT_NEAR(gt.r(0, 0), 0.32094, 1e-2) << "Expected first real coefficient to be close to 0.32094";
    ASSERT_NEAR(gt.i(0, 0), 0.0, 1e-2) << "Expected first imaginary coefficient to be close to 0.0";
    ASSERT_NEAR(gt.r(0, 1), 0.42759, 1e-2) << "Expected second real coefficient to be close to 0.42759";
    ASSERT_NEAR(gt.i(0, 1), 1.25295, 1e-2) << "Expected second imaginary coefficient to be close to 1.25295";
    ASSERT_NEAR(gt.r(0, 2), -0.44685, 1e-2) << "Expected third real coefficient to be close to -0.44685";
    ASSERT_NEAR(gt.i(0, 2), -0.65403, 1e-2) << "Expected third imaginary coefficient to be close to -0.65403";
}
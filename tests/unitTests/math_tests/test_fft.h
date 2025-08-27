#ifndef FFT_TEST_H
#define FFT_TEST_H

#include "../test_common.h"
#include <vector>
#include <cmath>

// Forward declarations for types used in FFT testing
struct fcType;

// Test fixture class for FFT functionality
class FFTTest : public ::testing::Test {
protected:
    // Set up test environment before each test
    void SetUp() override;
    
    // Clean up test environment after each test  
    void TearDown() override;
    
    // Helper methods for common test operations
    void CreateTemporalValues(int N, double x_start, double x_end, 
                             std::vector<std::vector<double>>& temporal_values);

    void InitializeFourierCoefficients(fcType& gt, int d = 1, int n = 16);
};

#endif // FFT_TEST_H
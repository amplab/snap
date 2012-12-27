#include "stdafx.h"
#include "Compat.h"
#include "TestLib.h"
#include "ProbabilityDistance.h"

// Test fixture for all the Landau-Viskhin Tests
struct ProbabilityDistanceTest {
    ProbabilityDistance dist;
    double prob;
    
    ProbabilityDistanceTest(): dist(0.1, 0.01, 0.2) {}
};

TEST_F(ProbabilityDistanceTest, "basic probabilities") {
    dist.compute("A", "A", "I", 1, 0, 0, &prob);
    ASSERT_NEAR(0.9, prob);

    dist.compute("A", "C", "I", 1, 0, 0, &prob);
    ASSERT_NEAR(0.1, prob);

    char quality10[2] = {43, 0};
    dist.compute("A", "C", quality10, 1, 0, 0, &prob);
    ASSERT_NEAR(0.19, prob);   // 1 - (1 - 0.9) * (1 - 0.9)

    dist.compute("A", "A", "I", 1, 1, 2, &prob);
    ASSERT_NEAR(0.9, prob);

    dist.compute("A", "C", "I", 1, 1, 2, &prob);
    ASSERT_NEAR(0.1, prob);

    dist.compute("A", "C", quality10, 1, 1, 2, &prob);
    ASSERT_NEAR(0.19, prob);   // 1 - (1 - 0.9) * (1 - 0.9)

    dist.compute("AAAAA", "AAAAA", "IIIII", 5, 1, 2, &prob);
    ASSERT_NEAR(pow(0.9, 5), prob);
}


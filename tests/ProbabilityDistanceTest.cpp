#include "stdafx.h"
#include "Compat.h"
#include "TestLib.h"
#include "ProbabilityDistance.h"

// Test fixture for all the ProbabilityDistance tests
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

    // Check that allowing a shift at the start doesn't change it
    dist.compute("A", "A", "I", 1, 1, 2, &prob);
    ASSERT_NEAR(0.9, prob);

    dist.compute("A", "C", "I", 1, 1, 2, &prob);
    ASSERT_NEAR(0.1, prob);

    dist.compute("A", "C", quality10, 1, 1, 2, &prob);
    ASSERT_NEAR(0.19, prob);   // 1 - (1 - 0.9) * (1 - 0.9)

    dist.compute("AAAAA", "AAAAA", "IIIII", 5, 1, 2, &prob);
    ASSERT_NEAR(pow(0.9, 5), prob);

    dist.compute("AAAAA", "AACAA", "IIIII", 5, 1, 2, &prob);
    ASSERT_NEAR(pow(0.9, 4) * 0.1, prob);
}


TEST_F(ProbabilityDistanceTest, "indels") {
    dist.compute("ACGTA", "ACGGTA", "IIIIII", 6, 1, 2, &prob);
    ASSERT_NEAR(pow(0.9, 5) * 0.01, prob);

    // Here it's better to count things as two substitutions than an indel and two mismatches
    dist.compute("ACGTA", "ACTA", "IIII", 4, 1, 2, &prob);
    ASSERT_NEAR(pow(0.9, 2) * pow(0.1, 2), prob);

    dist.compute("ACGTACGT", "ACGTTACGT", "IIIIIIIII", 9, 1, 2, &prob);
    ASSERT_NEAR(pow(0.9, 8) * 0.01, prob);

    dist.compute("ACGTACGT", "ACGACGT", "IIIIIII", 7, 1, 2, &prob);
    ASSERT_NEAR(pow(0.9, 7) * 0.01, prob);

    dist.compute("ACGTACGT", "ACTACGT", "IIIIIII", 7, 0, 2, &prob);
    ASSERT_NEAR(pow(0.9, 7) * 0.01, prob);

    // Here we can start at shift 1 and get a better probability with substitutions than indels
    dist.compute("ACGTACGT", "ACTACGT", "IIIIIII", 7, 1, 2, &prob);
    ASSERT_NEAR(pow(0.9, 5) * pow(0.1, 2), prob);

    dist.compute("ACGTACGT", "ACGTTTACGT", "IIIIIIIIII", 10, 1, 2, &prob);
    ASSERT_NEAR(pow(0.9, 8) * 0.01 * 0.2, prob);

    dist.compute("ACGTTTACGT", "ACGTACGT", "IIIIIIII", 8, 1, 2, &prob);
    ASSERT_NEAR(pow(0.9, 8) * 0.01 * 0.2, prob);
}

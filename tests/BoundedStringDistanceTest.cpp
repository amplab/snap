#include "stdafx.h"
#include "TestLib.h"
#include "BoundedStringDistance.h"

// Test fixture for all the Landau-Viskhin Tests
struct BoundedStringDistanceTest {
    BoundedStringDistance<> dist1, dist3, dist5, dist23, dist32;

    BoundedStringDistanceTest()
        : dist1(1, 1, 0.1, 0.01, 0.2), dist3(1, 3, 0.1, 0.01, 0.2), dist5(1, 5, 0.1, 0.01, 0.2),
          dist23(2, 3, 0.1, 0.01, 0.2), dist32(3, 2, 0.1, 0.01, 0.2)
    {}
};

TEST_F(BoundedStringDistanceTest, "equal strings") {
    ASSERT_EQ(0, dist1.compute("abcde", "abcde", 5, 0, 3));
    ASSERT_EQ(0, dist3.compute("abcde", "abcde", 5, 0, 3));
    ASSERT_EQ(0, dist5.compute("abcde", "abcde", 5, 0, 3));
    ASSERT_EQ(0, dist23.compute("abcde", "abcde", 5, 0, 3));
    ASSERT_EQ(0, dist32.compute("abcde", "abcde", 5, 0, 3));
}

TEST_F(BoundedStringDistanceTest, "prefixes") {
    ASSERT_EQ(0, dist1.compute("abcde", "abcd", 4, 0, 3));
    ASSERT_EQ(0, dist1.compute("abcde", "abc", 3, 0, 3));
    ASSERT_EQ(0, dist1.compute("abcde", "ab", 2, 0, 3));
    ASSERT_EQ(0, dist3.compute("abcde", "abcd", 4, 0, 3));
    ASSERT_EQ(0, dist3.compute("abcde", "abc", 3, 0, 3));
    ASSERT_EQ(0, dist3.compute("abcde", "ab", 2, 0, 3));
    ASSERT_EQ(0, dist5.compute("abcde", "abcd", 4, 0, 3));
    ASSERT_EQ(0, dist5.compute("abcde", "abc", 3, 0, 3));
    ASSERT_EQ(0, dist5.compute("abcde", "ab", 2, 0, 3));
    ASSERT_EQ(0, dist23.compute("abcde", "abcd", 4, 0, 3));
    ASSERT_EQ(0, dist23.compute("abcde", "abc", 3, 0, 3));
    ASSERT_EQ(0, dist23.compute("abcde", "ab", 2, 0, 3));
    ASSERT_EQ(0, dist32.compute("abcde", "abcd", 4, 0, 3));
    ASSERT_EQ(0, dist32.compute("abcde", "abc", 3, 0, 3));
    ASSERT_EQ(0, dist32.compute("abcde", "ab", 2, 0, 3));
}

TEST_F(BoundedStringDistanceTest, "gap open cost 1") {
    ASSERT_EQ(1, dist1.compute("abcdefg", "abcdXfg", 7, 0, 3));
    ASSERT_EQ(1, dist1.compute("abcdefg", "abdefg", 6, 0, 3));
    ASSERT_EQ(1, dist1.compute("abcdefg", "bcdefg", 6, 0, 3));
    ASSERT_EQ(2, dist1.compute("abcdefg", "cdefg", 5, 0, 3));
    ASSERT_EQ(3, dist1.compute("abcdefg", "defg", 4, 0, 3));
    ASSERT_EQ(1, dist1.compute("abcdefg", "abcXdefg", 8, 0, 3));
    ASSERT_EQ(2, dist1.compute("abcdefg", "abXXefg", 7, 0, 3));
    ASSERT_EQ(3, dist1.compute("abcdefg", "abXXXdefg", 9, 0, 3));
}

TEST_F(BoundedStringDistanceTest, "gap open cost 3") {
    ASSERT_EQ(1, dist3.compute("abcdefg", "abcdXfg", 7, 0, 3));
    ASSERT_EQ(3, dist3.compute("abcdefg", "abdefg", 6, 0, 3));
    ASSERT_EQ(3, dist3.compute("abcdefg", "bcdefg", 6, 0, 5));
    ASSERT_EQ(4, dist3.compute("abcdefg", "cdefg", 5, 0, 5));
    ASSERT_EQ(5, dist3.compute("abcdefghi", "aefghi", 6, 0, 5));
    ASSERT_EQ(3, dist3.compute("abcdefg", "abcXdefg", 8, 0, 5));
    ASSERT_EQ(4, dist3.compute("abcdefg", "abcXXdefg", 9, 0, 5));
    ASSERT_EQ(5, dist3.compute("abcdefg", "abcXXXdefg", 10, 0, 5));
}

TEST_F(BoundedStringDistanceTest, "sub cost 2, gap open cost 3") {
    ASSERT_EQ(2, dist23.compute("abcdefg", "abcdXfg", 7, 0, 3));
    ASSERT_EQ(4, dist23.compute("abcdefg", "XbcdXfg", 7, 0, 5));
    ASSERT_EQ(3, dist23.compute("abcdefg", "abdefg", 6, 0, 3));
    ASSERT_EQ(3, dist23.compute("abcdefg", "bcdefg", 6, 0, 5));
    ASSERT_EQ(4, dist23.compute("abcdefg", "cdefg", 5, 0, 5));
    ASSERT_EQ(5, dist23.compute("abcdefghi", "aefghi", 6, 0, 5));
    ASSERT_EQ(3, dist23.compute("abcdefg", "abcXdefg", 8, 0, 5));
    ASSERT_EQ(4, dist23.compute("abcdefg", "abcXXdefg", 9, 0, 5));
    ASSERT_EQ(5, dist23.compute("abcdefg", "abcXXXdefg", 10, 0, 5));
    // On this one it's better to count things as a deletion + insertion than 3 substitutions
    ASSERT_EQ(5, dist23.compute("abcdefg", "abcdXYZ", 7, 0, 6));
    ASSERT_EQ(6, dist23.compute("abcdefghi", "abcdXYZhi", 9, 0, 8));
}

TEST_F(BoundedStringDistanceTest, "sub cost 3, gap open cost 2") {
    ASSERT_EQ(3, dist32.compute("abcdefg", "abcdXfg", 7, 0, 3));
    ASSERT_EQ(6, dist32.compute("abcdefg", "XbcdXfg", 7, 0, 6));
    ASSERT_EQ(2, dist32.compute("abcdefg", "abdefg", 6, 0, 3));
    ASSERT_EQ(2, dist32.compute("abcdefg", "bcdefg", 6, 0, 5));
    ASSERT_EQ(3, dist32.compute("abcdefg", "cdefg", 5, 0, 5));
    ASSERT_EQ(4, dist32.compute("abcdefghi", "aefghi", 6, 0, 5));
    ASSERT_EQ(2, dist32.compute("abcdefg", "abcXdefg", 8, 0, 5));
    ASSERT_EQ(3, dist32.compute("abcdefg", "abcXXdefg", 9, 0, 5));
    ASSERT_EQ(4, dist32.compute("abcdefg", "abcXXXdefg", 10, 0, 5));
    // On this one it's better to count things as a deletion + insertion than 3 substitutions
    ASSERT_EQ(4, dist32.compute("abcdefg", "abcdXYZ", 7, 0, 6));
    ASSERT_EQ(8, dist32.compute("abcdefghijkl", "abcdXYZhijkl", 12, 0, 9));
}

TEST_F(BoundedStringDistanceTest, "overly distant strings") {
    ASSERT_EQ(-1, dist1.compute("abcde", "XXXXX", 5, 0, 2));
    ASSERT_EQ(-1, dist1.compute("abcdef", "def", 3, 0, 2));
    ASSERT_EQ(-1, dist1.compute("abcdef", "aXXXbcdef", 9, 0, 2));
    ASSERT_EQ(-1, dist3.compute("abcdef", "XXXXXX", 6, 0, 2));
    ASSERT_EQ(-1, dist3.compute("abcdefghi", "abdefghi", 8, 0, 2));
    ASSERT_EQ(-1, dist5.compute("abcdefghi", "abdefghi", 8, 0, 4));
    ASSERT_EQ(-1, dist23.compute("abcdefghi", "abdefghi", 8, 0, 2));
    ASSERT_EQ(-1, dist23.compute("abcdefghi", "abcdXfghi", 9, 0, 1));
    ASSERT_EQ(-1, dist32.compute("abcdefghi", "abdefghi", 8, 0, 1));
    ASSERT_EQ(-1, dist32.compute("abcdefghi", "abcdXfghi", 9, 0, 2));
}

TEST_F(BoundedStringDistanceTest, "CIGAR strings") {
    char cigarBuf[1024];
    int bufLen = sizeof(cigarBuf);

    ASSERT_EQ(0, dist3.compute("abcdefg", "abcdefg", 7, 0, 3, cigarBuf, bufLen));
    ASSERT_STREQ("7=", cigarBuf);

    ASSERT_EQ(1, dist3.compute("abcdefg", "abcdXfg", 7, 0, 3, cigarBuf, bufLen));
    ASSERT_STREQ("4=1X2=", cigarBuf);

    ASSERT_EQ(1, dist3.compute("abcdefg", "Xbcdefg", 7, 0, 3, cigarBuf, bufLen));
    ASSERT_STREQ("1X6=", cigarBuf);

    ASSERT_EQ(1, dist3.compute("abcdefg", "abcdefG", 7, 0, 3, cigarBuf, bufLen));
    ASSERT_STREQ("6=1X", cigarBuf);

    ASSERT_EQ(3, dist3.compute("abcdefg", "abdefg", 6, 0, 3, cigarBuf, bufLen));
    ASSERT_STREQ("2=1D4=", cigarBuf);

    ASSERT_EQ(3, dist3.compute("abcdefg", "bcdefg", 6, 0, 5, cigarBuf, bufLen));
    ASSERT_STREQ("1D6=", cigarBuf);

    ASSERT_EQ(4, dist3.compute("abcdefg", "cdefg", 5, 0, 5, cigarBuf, bufLen));
    ASSERT_STREQ("2D5=", cigarBuf);

    ASSERT_EQ(5, dist3.compute("abcdefghi", "aefghi", 6, 0, 5, cigarBuf, bufLen));
    ASSERT_STREQ("1=5X", cigarBuf);

    ASSERT_EQ(3, dist3.compute("abcdefg", "abcXdefg", 8, 0, 5, cigarBuf, bufLen));
    ASSERT_STREQ("3=1I4=", cigarBuf);

    ASSERT_EQ(4, dist3.compute("abcdefg", "abcXXdefg", 9, 0, 5, cigarBuf, bufLen));
    ASSERT_STREQ("3=2I4=", cigarBuf);

    ASSERT_EQ(5, dist3.compute("abcdefg", "abcXXXdefg", 10, 0, 5, cigarBuf, bufLen));
    ASSERT_STREQ("3=3I4=", cigarBuf);

    ASSERT_EQ(6, dist3.compute("abcdefghij", "bcdefgXhij", 10, 0, 6, cigarBuf, bufLen));
    ASSERT_STREQ("1D6=1I3=", cigarBuf);

    // Versions with useM = true

    ASSERT_EQ(0, dist3.compute("abcdefg", "abcdefg", 7, 0, 3, cigarBuf, bufLen, true));
    ASSERT_STREQ("7M", cigarBuf);

    ASSERT_EQ(1, dist3.compute("abcdefg", "abcdXfg", 7, 0, 3, cigarBuf, bufLen, true));
    ASSERT_STREQ("7M", cigarBuf);

    ASSERT_EQ(6, dist3.compute("abcdefghij", "bcdefgXhij", 10, 0, 6, cigarBuf, bufLen, true));
    ASSERT_STREQ("1D6M1I3M", cigarBuf);

    // Versions with substitutionPenalty > 1
    
    ASSERT_EQ(4, dist32.compute("abcdefg", "abcdXYZ", 7, 0, 6, cigarBuf, bufLen));
    ASSERT_STREQ("4=3I", cigarBuf);

    ASSERT_EQ(6, dist32.compute("abcdefghi", "abcdXYZhi", 9, 0, 6, cigarBuf, bufLen));
    ASSERT_STREQ("4=5I", cigarBuf);

    ASSERT_EQ(8, dist32.compute("abcdefghijkl", "abcdXYZhijkl", 12, 0, 9, cigarBuf, bufLen));
    ASSERT_STREQ("4=3I3D5=", cigarBuf);
}

TEST_F(BoundedStringDistanceTest, "map probability") {
    double prob;
    const char *highQuality = "IIIIIIIIIIIIIIIIII";
    const char quality10[2] = {43, 0};

    ASSERT_EQ(0, dist23.compute("a", "a", 1, 0, 3, &prob, highQuality));
    ASSERT_NEAR(0.9, prob);

    ASSERT_EQ(2, dist23.compute("a", "X", 1, 0, 3, &prob, highQuality));
    ASSERT_NEAR(0.1, prob);

    ASSERT_EQ(2, dist23.compute("a", "X", 1, 0, 3, &prob, quality10));
    ASSERT_NEAR(0.19, prob);    // 1 - (1 - 0.9) * (1 - 0.9)

    ASSERT_EQ(0, dist23.compute("aaaaa", "aaaaa", 5, 0, 3, &prob, highQuality));
    ASSERT_NEAR(pow(0.9, 5), prob);

    ASSERT_EQ(2, dist23.compute("aaaaa", "aaXaa", 5, 0, 3, &prob, highQuality));
    ASSERT_NEAR(pow(0.9, 4) * 0.1, prob);

    ASSERT_EQ(3, dist23.compute("aaaaa", "aaXaaa", 6, 0, 3, &prob, highQuality));
    ASSERT_NEAR(pow(0.9, 5) * 0.01, prob);

    ASSERT_EQ(4, dist23.compute("aaaaa", "aaXXaaa", 7, 0, 5, &prob, highQuality));
    ASSERT_NEAR(pow(0.9, 5) * 0.01 * 0.2, prob);

    ASSERT_EQ(3, dist23.compute("abcdef", "acdef", 5, 0, 3, &prob, highQuality));
    ASSERT_NEAR(pow(0.9, 5) * 0.01, prob);

    ASSERT_EQ(4, dist23.compute("abcdefg", "adefg", 5, 0, 5, &prob, highQuality));
    ASSERT_NEAR(pow(0.9, 5) * 0.01 * 0.2, prob);

    ASSERT_EQ(4, dist23.compute("abcdefg", "cdefg", 5, 0, 5, &prob, highQuality));
    ASSERT_NEAR(pow(0.9, 5) * 0.01 * 0.2, prob);

    ASSERT_EQ(5, dist23.compute("abcdefg", "abcdXYZ", 7, 0, 6, &prob, highQuality));
    ASSERT_NEAR(pow(0.9, 4) * 0.01 * 0.2 * 0.2, prob);  // This will be counted as one long indel

    ASSERT_EQ(7, dist23.compute("abcdefg", "XbcdXYZ", 7, 0, 7, &prob, highQuality));
    ASSERT_NEAR(pow(0.9, 3) * 0.1 * 0.01 * 0.2 * 0.2, prob);  // Indel plus substitution

    ASSERT_EQ(6, dist23.compute("abcdefghi", "abcdXYZhi", 9, 0, 8, &prob, highQuality));
    ASSERT_NEAR(pow(0.9, 6) * pow(0.1, 3), prob);
}

TEST_F(BoundedStringDistanceTest, "start shifts") {
    ASSERT_EQ(3, dist3.compute("abcdefg", "bcdefg", 6, 0, 5));
    ASSERT_EQ(0, dist3.compute("abcdefg", "bcdefg", 6, 1, 5));
    ASSERT_EQ(0, dist3.compute("abcdefg", "bcdefg", 6, 2, 5));
    ASSERT_EQ(4, dist3.compute("abcdefg", "bcdXfg", 6, 0, 5));
    ASSERT_EQ(1, dist3.compute("abcdefg", "bcdXfg", 6, 1, 5));
    ASSERT_EQ(1, dist3.compute("abcdefg", "bcdXfg", 6, 2, 5));

    ASSERT_EQ(3, dist3.compute("abcdefg" + 2, "bcdefg", 6, 0, 5));
    ASSERT_EQ(0, dist3.compute("abcdefg" + 2, "bcdefg", 6, 1, 5));
    ASSERT_EQ(0, dist3.compute("abcdefg" + 2, "bcdefg", 6, 2, 5));
    ASSERT_EQ(4, dist3.compute("abcdefg" + 2, "bcdXfg", 6, 0, 5));
    ASSERT_EQ(1, dist3.compute("abcdefg" + 2, "bcdXfg", 6, 1, 5));
    ASSERT_EQ(1, dist3.compute("abcdefg" + 2, "bcdXfg", 6, 2, 5));

    ASSERT_EQ(4, dist3.compute("abcdefg", "cdefg", 5, 0, 5));
    ASSERT_EQ(3, dist3.compute("abcdefg", "cdefg", 5, 1, 5));
    ASSERT_EQ(0, dist3.compute("abcdefg", "cdefg", 5, 2, 5));
}

TEST_F(BoundedStringDistanceTest, "empty strings") {
    double prob;
    ASSERT_EQ(0, dist23.compute("a", "", 0, 0, 3, &prob, "I"));
    ASSERT_NEAR(1.0, prob);

    char cigarBuf[1024];
    int bufLen = sizeof(cigarBuf);
    ASSERT_EQ(0, dist23.compute("abcdefg", "", 0, 0, 3, cigarBuf, bufLen));
    ASSERT_STREQ("", cigarBuf);
}

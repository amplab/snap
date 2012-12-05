#include "stdafx.h"
#include "TestLib.h"
#include "BoundedStringDistance.h"

// Test fixture for all the Landau-Viskhin Tests
struct BoundedStringDistanceTest {
    BoundedStringDistance<> dist1, dist3, dist5;
    BoundedStringDistance<true> dist1Cigar, dist3Cigar;

    BoundedStringDistanceTest()
        : dist1(1), dist3(3), dist5(5), dist1Cigar(1), dist3Cigar(3)
    {}
};

TEST_F(BoundedStringDistanceTest, "equal strings") {
    ASSERT_EQ(0, dist1.compute("abcde", "abcde", 5, 3));
    ASSERT_EQ(0, dist3.compute("abcde", "abcde", 5, 3));
    ASSERT_EQ(0, dist5.compute("abcde", "abcde", 5, 3));
}

TEST_F(BoundedStringDistanceTest, "prefixes") {
    ASSERT_EQ(0, dist1.compute("abcde", "abcd", 4, 3));
    ASSERT_EQ(0, dist1.compute("abcde", "abc", 3, 3));
    ASSERT_EQ(0, dist1.compute("abcde", "ab", 2, 3));
    ASSERT_EQ(0, dist3.compute("abcde", "abcd", 4, 3));
    ASSERT_EQ(0, dist3.compute("abcde", "abc", 3, 3));
    ASSERT_EQ(0, dist3.compute("abcde", "ab", 2, 3));
    ASSERT_EQ(0, dist5.compute("abcde", "abcd", 4, 3));
    ASSERT_EQ(0, dist5.compute("abcde", "abc", 3, 3));
    ASSERT_EQ(0, dist5.compute("abcde", "ab", 2, 3));
}

TEST_F(BoundedStringDistanceTest, "gap open cost 1") {
    ASSERT_EQ(1, dist1.compute("abcdefg", "abcdXfg", 7, 3));
    ASSERT_EQ(1, dist1.compute("abcdefg", "abdefg", 6, 3));
    ASSERT_EQ(1, dist1.compute("abcdefg", "bcdefg", 6, 3));
    ASSERT_EQ(2, dist1.compute("abcdefg", "cdefg", 5, 3));
    ASSERT_EQ(3, dist1.compute("abcdefg", "defg", 4, 3));
    ASSERT_EQ(1, dist1.compute("abcdefg", "abcXdefg", 8, 3));
    ASSERT_EQ(2, dist1.compute("abcdefg", "abXXefg", 7, 3));
    ASSERT_EQ(3, dist1.compute("abcdefg", "abXXXdefg", 9, 3));
}

TEST_F(BoundedStringDistanceTest, "gap open cost 3") {
    ASSERT_EQ(1, dist3.compute("abcdefg", "abcdXfg", 7, 3));
    ASSERT_EQ(3, dist3.compute("abcdefg", "abdefg", 6, 3));
    ASSERT_EQ(3, dist3.compute("abcdefg", "bcdefg", 6, 5));
    ASSERT_EQ(4, dist3.compute("abcdefg", "cdefg", 5, 5));
    ASSERT_EQ(5, dist3.compute("abcdefghi", "aefghi", 6, 5));
    ASSERT_EQ(3, dist3.compute("abcdefg", "abcXdefg", 8, 5));
    ASSERT_EQ(4, dist3.compute("abcdefg", "abcXXdefg", 9, 5));
    ASSERT_EQ(5, dist3.compute("abcdefg", "abcXXXdefg", 10, 5));
}

TEST_F(BoundedStringDistanceTest, "overly distant strings") {
    ASSERT_EQ(-1, dist1.compute("abcde", "XXXXX", 5, 2));
    ASSERT_EQ(-1, dist1.compute("abcdef", "def", 3, 2));
    ASSERT_EQ(-1, dist1.compute("abcdef", "aXXXbcdef", 9, 2));
    ASSERT_EQ(-1, dist3.compute("abcdef", "XXXXXX", 6, 2));
    ASSERT_EQ(-1, dist3.compute("abcdefghi", "abdefghi", 8, 2));
    ASSERT_EQ(-1, dist5.compute("abcdefghi", "abdefghi", 8, 4));
}

TEST_F(BoundedStringDistanceTest, "CIGAR strings") {
    char cigarBuf[1024];
    int bufLen = sizeof(cigarBuf);

    ASSERT_EQ(0, dist3Cigar.compute("abcdefg", "abcdefg", 7, 3, cigarBuf, bufLen));
    ASSERT_STREQ("7=", cigarBuf);

    ASSERT_EQ(1, dist3Cigar.compute("abcdefg", "abcdXfg", 7, 3, cigarBuf, bufLen));
    ASSERT_STREQ("4=1X2=", cigarBuf);

    ASSERT_EQ(1, dist3Cigar.compute("abcdefg", "Xbcdefg", 7, 3, cigarBuf, bufLen));
    ASSERT_STREQ("1X6=", cigarBuf);

    ASSERT_EQ(1, dist3Cigar.compute("abcdefg", "abcdefG", 7, 3, cigarBuf, bufLen));
    ASSERT_STREQ("6=1X", cigarBuf);

    ASSERT_EQ(3, dist3Cigar.compute("abcdefg", "abdefg", 6, 3, cigarBuf, bufLen));
    ASSERT_STREQ("2=1D4=", cigarBuf);

    ASSERT_EQ(3, dist3Cigar.compute("abcdefg", "bcdefg", 6, 5, cigarBuf, bufLen));
    ASSERT_STREQ("1D6=", cigarBuf);

    ASSERT_EQ(4, dist3Cigar.compute("abcdefg", "cdefg", 5, 5, cigarBuf, bufLen));
    ASSERT_STREQ("2D5=", cigarBuf);

    ASSERT_EQ(5, dist3Cigar.compute("abcdefghi", "aefghi", 6, 5, cigarBuf, bufLen));
    ASSERT_STREQ("1=5X", cigarBuf);

    ASSERT_EQ(3, dist3Cigar.compute("abcdefg", "abcXdefg", 8, 5, cigarBuf, bufLen));
    ASSERT_STREQ("3=1I4=", cigarBuf);

    ASSERT_EQ(4, dist3Cigar.compute("abcdefg", "abcXXdefg", 9, 5, cigarBuf, bufLen));
    ASSERT_STREQ("3=2I4=", cigarBuf);

    ASSERT_EQ(5, dist3Cigar.compute("abcdefg", "abcXXXdefg", 10, 5, cigarBuf, bufLen));
    ASSERT_STREQ("3=3I4=", cigarBuf);

    ASSERT_EQ(6, dist3Cigar.compute("abcdefghij", "bcdefgXhij", 10, 6, cigarBuf, bufLen));
    ASSERT_STREQ("1D6=1I3=", cigarBuf);

    // Versions with useM = true

    ASSERT_EQ(0, dist3Cigar.compute("abcdefg", "abcdefg", 7, 3, cigarBuf, bufLen, true));
    ASSERT_STREQ("7M", cigarBuf);

    ASSERT_EQ(1, dist3Cigar.compute("abcdefg", "abcdXfg", 7, 3, cigarBuf, bufLen, true));
    ASSERT_STREQ("7M", cigarBuf);

    ASSERT_EQ(6, dist3Cigar.compute("abcdefghij", "bcdefgXhij", 10, 6, cigarBuf, bufLen, true));
    ASSERT_STREQ("1D6M1I3M", cigarBuf);
}

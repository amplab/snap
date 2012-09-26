#include "stdafx.h"
#include "TestLib.h"
#include "LandauVishkin.h"

// Test fixture for all the Landau-Viskhin Tests
struct LandauVishkinTest {
    LandauVishkin lv;
    LandauVishkinWithCigar lvc;
};

TEST_F(LandauVishkinTest, "equal strings") {
    ASSERT_EQ(0, lv.computeEditDistance("abcde", 5, "abcde", 5, 2));
}

TEST_F(LandauVishkinTest, "prefixes") {
    ASSERT_EQ(0, lv.computeEditDistance("abcde", 5, "abcd", 4, 2));
    ASSERT_EQ(0, lv.computeEditDistance("abcde", 5, "abc", 3, 2));
    ASSERT_EQ(0, lv.computeEditDistance("abcde", 5, "ab", 2, 2));
}

TEST_F(LandauVishkinTest, "non-equal strings") {
    ASSERT_EQ(1, lv.computeEditDistance("abcde", 5, "abcdX", 5, 2));
    ASSERT_EQ(1, lv.computeEditDistance("abcde", 5, "abde", 4, 2));
    ASSERT_EQ(1, lv.computeEditDistance("abcde", 5, "bcde", 4, 2));
    ASSERT_EQ(1, lv.computeEditDistance("abcde", 5, "abcXde", 6, 2));
    ASSERT_EQ(2, lv.computeEditDistance("abcde", 5, "abXXe", 5, 2));
    ASSERT_EQ(2, lv.computeEditDistance("abcde", 5, "abcXXde", 7, 2));
}

TEST_F(LandauVishkinTest, "overly distant strings") {
    ASSERT_EQ(-1, lv.computeEditDistance("abcde", 5, "XXXXX", 5, 2));
}

TEST_F(LandauVishkinTest, "CIGAR strings") {
    char cigarBuf[1024];
    int bufLen = sizeof(cigarBuf);

    lvc.computeEditDistance("abcde", 5, "abcde", 5, 2, cigarBuf, bufLen, false);
    ASSERT_STREQ("5=", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "abcde", 5, 2, cigarBuf, bufLen, true);
    ASSERT_STREQ("5M", cigarBuf);

    lvc.computeEditDistance("abcdef", 6, "abcde", 5, 2, cigarBuf, bufLen, false);
    ASSERT_STREQ("5=", cigarBuf);

    lvc.computeEditDistance("abcdef", 6, "abcde", 5, 2, cigarBuf, bufLen, true);
    ASSERT_STREQ("5M", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "abcdX", 5, 2, cigarBuf, bufLen, false);
    ASSERT_STREQ("4=1X", cigarBuf);   // This used to give 4=1I before

    lvc.computeEditDistance("abcde", 5, "abcdX", 5, 2, cigarBuf, bufLen, true);
    ASSERT_STREQ("5M", cigarBuf);   // This used to give 4=1I before

    lvc.computeEditDistance("abcde", 5, "Xbcde", 5, 2, cigarBuf, bufLen, false);
    ASSERT_STREQ("1X4=", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "Xbcde", 5, 2, cigarBuf, bufLen, true);
    ASSERT_STREQ("5M", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "abde", 4, 2, cigarBuf, bufLen, false);
    ASSERT_STREQ("2=1D2=", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "abde", 4, 2, cigarBuf, bufLen, true);
    ASSERT_STREQ("2M1D2M", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "bcde", 4, 2, cigarBuf, bufLen, false);
    ASSERT_STREQ("1D4=", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "bcde", 4, 2, cigarBuf, bufLen, true);
    ASSERT_STREQ("1D4M", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "abcXde", 6, 2, cigarBuf, bufLen, false);
    ASSERT_STREQ("3=1I2=", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "abcXde", 6, 2, cigarBuf, bufLen, true);
    ASSERT_STREQ("3M1I2M", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "abXXe", 5, 2, cigarBuf, bufLen, false);
    ASSERT_STREQ("2=2X1=", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "abXXe", 5, 2, cigarBuf, bufLen, true);
    ASSERT_STREQ("5M", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "abcXXde", 7, 3, cigarBuf, bufLen, false);
    ASSERT_STREQ("3=2I2=", cigarBuf);

    lvc.computeEditDistance("abcde", 5, "abcXXde", 7, 3, cigarBuf, bufLen, true);
    ASSERT_STREQ("3M2I2M", cigarBuf);

    lvc.computeEditDistance("ttttc", 5, "tttc", 4, 3, cigarBuf, bufLen, false);
    ASSERT_STREQ("3=1X", cigarBuf);

    lvc.computeEditDistance("ttttc", 5, "tttc", 4, 3, cigarBuf, bufLen, true);
    ASSERT_STREQ("4M", cigarBuf);

    lvc.computeEditDistance("tttcc", 5, "ttttc", 5, 3, cigarBuf, bufLen, false);
    ASSERT_STREQ("3=1X1=", cigarBuf);

    lvc.computeEditDistance("tttcc", 5, "ttttc", 5, 3, cigarBuf, bufLen, true);
    ASSERT_STREQ("5M", cigarBuf);

    lvc.computeEditDistance("tttcc", 5, "tttaa", 5, 3, cigarBuf, bufLen, false);
    ASSERT_STREQ("3=2X", cigarBuf);

    lvc.computeEditDistance("tttcc", 5, "tttaa", 5, 3, cigarBuf, bufLen, true);
    ASSERT_STREQ("5M", cigarBuf);

    // A real example where we used to give 1D and later 1I instead of 2X
    lvc.computeEditDistance("atctcag", 7, "acttcag", 7, 3, cigarBuf, bufLen, false);
    ASSERT_STREQ("1=2X4=", cigarBuf);

    lvc.computeEditDistance("atctcag", 7, "acttcag", 7, 3, cigarBuf, bufLen, true);
    ASSERT_STREQ("7M", cigarBuf);

    // Edge cases when pattern is longer than available text
    lvc.computeEditDistance("abc", 3, "abcde", 5, 3, cigarBuf, bufLen, false);
    ASSERT_STREQ("3=2X", cigarBuf);

    lvc.computeEditDistance("abc", 3, "abcde", 5, 3, cigarBuf, bufLen, true);
    ASSERT_STREQ("5M", cigarBuf);

    lvc.computeEditDistance("abc", 3, "abXde", 5, 3, cigarBuf, bufLen, false);
    ASSERT_STREQ("2=3X", cigarBuf);

    lvc.computeEditDistance("abc", 3, "abXde", 5, 3, cigarBuf, bufLen, true);
    ASSERT_STREQ("5M", cigarBuf);
}

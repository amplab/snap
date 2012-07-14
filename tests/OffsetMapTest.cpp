#include "../SNAPLib/OffsetMap.h"
#include "TestLib.h"

struct OffsetMapTest {
    void addMutation(unsigned loc, unsigned before, unsigned after) {
        offsetMap.addMutation(loc, before, after);
    }
    void test(unsigned expectedAncestor, unsigned descendant) {
        ASSERT_EQ(expectedAncestor, offsetMap.getTranslation(descendant));
    }
    sim::OffsetMap offsetMap;
};

TEST_F(OffsetMapTest, "trivial") {
    for (unsigned i = 0; i < 100; ++i)
        test(i, i);
}

TEST_F(OffsetMapTest, "delete 3") {
    addMutation(0, 3, 0);
    test(3, 0);
    test(4, 1);
}

TEST_F(OffsetMapTest, "insert 1") {
    addMutation(0, 1, 2);
    test(0, 0);
    test(0, 1);
    test(1, 2);
}

TEST_F(OffsetMapTest, "insert 2") {
    addMutation(1, 1, 3);
    test(0, 0);
    test(1, 1);
    test(1, 2);
    test(1, 3);
    test(2, 4);
}

TEST_F(OffsetMapTest, "insert and delete") {
    addMutation(1, 2, 3);
    addMutation(6, 3, 0);
    addMutation(10, 2, 1);
    test(2, 3);
    test(3, 4);
    test(4, 5);
    test(10, 8);
    test(12, 9);
}

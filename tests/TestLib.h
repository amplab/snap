#pragma once

/**
 * A tiny unit testing library in the spirit of Google Test.
 *
 * To inplement a standalone test, write:
 *
 *    TEST("description") { body }
 *
 * For fixture-based tests, define a struct Fixture with the fields you want
 * available (all public) and any setup and teardown code, then use TEST_F:
 *
 *    struct MyFixture {
 *        int field1;
 *        MyFixture() {}  // Optional setup code
 *        ~MyFixture() {} // Optional teardown code
 *    }
 *    
 *    TEST_F(MyFixture, "description") { body }
 *
 *    TEST_F(MyFixture, "description 2") { another body }
 *
 * In the body of a test, you can use the following macros and assertions:
 *
 *    ASSERT(expression)
 *    ASSERT_M(expression, message)
 *    ASSERT_EQ(expected, actualValue)
 *    ASSERT_NE(expected, actualValue)
 *    ASSERT_STREQ(expected, actualValue)   (for C strings)
 *    ASSERT_STRNE(expected, actualValue)
 *    FAIL(message)
 */

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

namespace test {

class TestCase;

typedef void (*FunctionPtr)();

struct TestCase {
    TestCase(const char *fixture_, const char *name_, FunctionPtr func_)
            : fixture(fixture_), name(name_), func(func_) {
        getCases().push_back(this);
    }
    
    void run() { func(); };

    const char *fixture;
    const char *name;
    FunctionPtr func;

    static std::vector<TestCase*>& getCases() {
        static std::vector<TestCase*> cases;
        return cases;
    };
};

struct TestFailedException {
    TestFailedException(const char *file_, int line_, const std::string& message_)
        : file(file_), line(line_), message(message_) {}

    const char *file;
    int line;
    std::string message;
};

int runAllTests();

}

#define CONCAT1( x, y ) x ## y
#define CONCAT2( x, y ) CONCAT1( x, y ) /* To escape weird macro expansion rules */
#define TEST_FUNC(line)  CONCAT2(_test_func_,  line)
#define TEST_CASE(line)  CONCAT2(_test_case_,  line)
#define TEST_CLASS(line) CONCAT2(_test_class_, line)

#define TEST(name) \
    static void TEST_FUNC(__LINE__) (); \
    static test::TestCase TEST_CASE(__LINE__) (__FILE__, name, &TEST_FUNC(__LINE__)); \
    static void TEST_FUNC(__LINE__) () /* body follows */

#define TEST_F(fixture, name) \
    namespace { struct TEST_CLASS(__LINE__) : public fixture { void _run(); }; } \
    static void TEST_FUNC(__LINE__) () { TEST_CLASS(__LINE__) cls; cls._run(); } \
    static test::TestCase TEST_CASE(__LINE__) (#fixture, name, &TEST_FUNC(__LINE__)); \
    void TEST_CLASS(__LINE__)::_run() /* body follows */

#define ASSERT(expr) \
    if (!(expr)) { \
        std::ostringstream oss; \
        oss << "assertion failed: " #expr; \
        throw test::TestFailedException(__FILE__, __LINE__, oss.str()); \
    }

#define ASSERT_M(expr, message) \
    if (!(expr)) { \
        std::ostringstream oss; \
        oss << "assertion failed: " << message; \
        throw test::TestFailedException(__FILE__, __LINE__, oss.str()); \
    }

#define ASSERT_EQ(expected, actual) \
    if (!((expected) == (actual))) { \
        std::ostringstream oss; \
        oss << #actual << " was " << (actual) << ", expected " << (expected); \
        throw test::TestFailedException(__FILE__, __LINE__, oss.str()); \
    }

#define ASSERT_NE(expected, actual) \
    if (!((expected) != (actual))) { \
        std::ostringstream oss; \
        oss << #actual << " was " << (expected); \
        throw test::TestFailedException(__FILE__, __LINE__, oss.str()); \
    }

#define ASSERT_STREQ(expected, actual) \
    if (strcmp((expected), (actual)) != 0) { \
        std::ostringstream oss; \
        oss << #actual << " was \"" << (actual) << "\", expected \"" << (expected) << "\""; \
        throw test::TestFailedException(__FILE__, __LINE__, oss.str()); \
    }

#define ASSERT_STRNE(expected, actual) \
    if (strcmp((expected), (actual)) == 0) { \
        std::ostringstream oss; \
        oss << #actual << " was \"" << (expected) << "\""; \
        throw test::TestFailedException(__FILE__, __LINE__, oss.str()); \
    }

#define FAIL(message) \
    throw test::TestFailedException(__FILE__, __LINE__, message);

#include <iostream>
#include <cstring>

#include "TestLib.h"

using namespace std;
using namespace test;

int test::runAllTests() {
    const std::vector<TestCase*> &testCases = TestCase::getCases();
    int passed = 0;
    const char *prevFixture = "";

    for (int i = 0; i < testCases.size(); i++) {
        TestCase *tc = testCases[i];
        if (strcmp(tc->fixture, prevFixture) != 0) {
            if (strlen(prevFixture) != 0) {
                cout << endl;
            }
            cout << tc->fixture << ":" << endl;
            prevFixture = tc->fixture;
        }
        cout << "- " << tc->name << ": " << flush;
        try {
            tc->run();
            cout << "[OK]" << endl;
            passed++;
        } catch (TestFailedException &e) {
            cout << "[FAILED]" << endl;
            cout << "    " << e.message << endl;
            cout << "    (" << e.file << ":" << e.line << ")" << endl;
        }
    }
    
    cout << endl << passed << " / " << testCases.size() << " tests passed." << endl;
    return (passed == testCases.size() ? 0 : 1);
}

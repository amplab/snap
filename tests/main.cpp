#include "TestLib.h"

int main(int argc, char **argv) {
    // Allow passing in a substring to search for in test names
    char *filter = (argc == 2 ? argv[1] : NULL);
    return test::runAllTests(filter);
}

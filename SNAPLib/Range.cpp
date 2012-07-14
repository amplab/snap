#include "stdafx.h"
#include "Range.h"

using namespace std;


Range Range::parse(const char *str)
{
    size_t len = strlen(str);
    int colons = 0;
    for (int i = 0; i < len; i++) {
        if (str[i] == ':') {
            colons++;
        } else if (!isdigit(str[i])) {
            fprintf(stderr, "Invalid format: %s\n", str);
            exit(1);
        }
    }
    int start, end, step;
    if (colons == 0) {
        sscanf(str, "%d", &start);
        end = start;
        step = 1;
    } else if (colons == 1) {
        sscanf(str, "%d:%d", &start, &end);
        step = 1;
    } else if (colons == 2) {
        sscanf(str, "%d:%d:%d", &start, &step, &end);
    } else {
        fprintf(stderr, "Invalid format: %s\n", str);
        exit(1);
    }
    if (step < 1 || end < start) {
        fprintf(stderr, "Invalid range: %s\n", str);
        exit(1);
    }
    return Range(start, end, step);
}

#pragma once


//
// An inclusive range class, used mostly to iterate over algorithm parameter settings.
//
class Range
{
public:
    int start, end, step;

    Range() : start(0), end(0), step(1) {}
    Range(int number) : start(number), end(number), step(1) {}
    Range(int start_, int end_) : start(start_), end(end_), step(1) {}
    Range(int start_, int end_, int step_) : start(start_), end(end_), step(step_) {}

    int size() { return 1 + (end - start) / step; }

    // 
    // Parses a range passed as a string. Supports the following formats:
    // 
    // num       single number
    // n1:n2     numbers from n1 to n2 inclusive, stepping by 1
    // n1:s:n2   numbers from n1 to n2 inclusive, stepping by s
    //
    static Range parse(const char *str);
};


//
// Macro for iterating over a Range.
//
#define FOR_RANGE(elem, range) for (int elem = range.start; elem <= range.end; elem += range.step)

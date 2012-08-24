#pragma once

#include <map>
#include "stdafx.h"
#include "Compat.h"
#include "GoodRandom.h"
#include "Tables.h"

//
// General utilities.
//

namespace util {

#define BEGEND(container) (container).begin(), (container).end()

//
// Sample from the Bernoulli distribution corresponding to 'probability'.
// 
inline bool sampleBernoulli(double probability, _uint64 resolution=0xffffffff)
{
    _ASSERT(0 <= probability && probability <= 1);
    return GoodFastRandom(resolution) / (double)resolution < probability;
}

inline double ratio(double a, double b=1)
{
    return a / (a + b);
}

inline void die(const std::string &epitaph, int code=1)
{
    fputs((epitaph + '\n').c_str(), stderr);
    exit(code);
}

const int MAXLINE = 1024;

inline bool readLine(char *line, FILE *file)
{
    return fgets(line, MAXLINE, file) && strchr(line, '\n');
}

inline void dieAtLine(const std::string &epitaph, const std::string &line)
{
    die(epitaph + "\n  line: " + line);
}

//
// Read a line, returning whether the read was successful or EOF occurred.
//
// Die if an error occurred.
// 
inline bool readLineOrDie(char *line, FILE *file)
{
    if (!readLine(line, file)) {
        if (feof(file))
            return false;
        dieAtLine("Error reading line.", line);
    }
    return true;
}

template <class Key, class T>
const T &getOrElse(const std::map<Key, T> &m, typename std::map<Key, T>::const_iterator p,
                   const T &d=T())
{
    return p == m.end() ? d : p->second;
}

//
// Scala map method: lookup with default if not found.
// 
template <class Key, class T>
const T &getOrElse(const std::map<Key, T> &m, const Key &k, const T &d=T())
{
    return getOrElse(m, m.find(k), d);
}

template <class RandIter, class T>
size_t findIndex(RandIter first, RandIter last, const T &value)
{
    return std::find(first, last, value) - first;
}

//
// Analogue of Scala's 'mkString' method.
// 
template <class InIter>
std::string joinWithSep(InIter first, InIter last, char sep)
{
    std::string joined;
    if (first == last)
        return joined;
    joined += *first++;
    while (first != last)
        (joined += sep) += *first++;
    return joined;
}

template <class T>
void addIfAbsent(std::vector<T> *v, const T &e)
{
    if (!std::count(v->begin(), v->end(), e))
        v->push_back(e);
}

inline FILE *fopenOrDie(const char *fileName, const char *mode)
{
    FILE *file = fopen(fileName, mode);
    if (!file)
        die("Unable to open file '%s'.");
    return file;
}

inline void fcloseOrDie(FILE *file)
{
    if (ferror(file))
        die("Error in file about to be closed.");
    fclose(file);
}

inline bool startsWith(const char *s, const char *prefix)
{
    return strstr(s, prefix) == s;
}

inline const char *strstrAfter(const char *s, const char *t)
{
    const char *tins = strstr(s, t);
    return tins ? tins + strlen(t) : NULL;
}

    inline void
toComplement(
    char* rc,
    const char* bases = NULL,
    int length = -1,
    bool toLower = false)
{
    if (length < 0) {
        length = (int) strlen(bases != NULL ? bases : rc);
        if (bases != NULL) {
            rc[length] = '\0';
        }
    }
    if (bases != NULL) {
        if (! toLower) {
            for (int i = 0; i < length; i++) {
                rc[i] = COMPLEMENT[bases[length - i - 1]];
            }
        } else {
            for (int i = 0; i < length; i++) {
                rc[i] = tolower(COMPLEMENT[bases[length - i - 1]]);
            }
        }
    } else {
        // reverse complement in place
        for (int i = 0; i < length / 2; i++) {
            char t = COMPLEMENT[rc[i]];
            rc[i] = toLower ? tolower(COMPLEMENT[rc[length - i - 1]]) : COMPLEMENT[rc[length - i - 1]];
            rc[length - i - 1] = toLower ? tolower(t) : t;
        }
        if (length % 2 == 1) {
            rc[length / 2] = toLower ? tolower(COMPLEMENT[rc[length / 2]]) : COMPLEMENT[rc[length / 2]];
        }
    }
}

    
    inline unsigned
log10bucket(
    unsigned n)
{
    unsigned factor = 1;
    while (n >= 10) {
        n /= 10;
        factor *= 10;
    }
    return n * factor;
}

} // namespace util

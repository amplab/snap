#pragma once

#include <map>
#include "stdafx.h"
#include "Compat.h"
#include "Tables.h"
#include "exit.h"

//
// General utilities.
//

namespace util {

#define BEGEND(container) (container).begin(), (container).end()

inline double ratio(double a, double b=1)
{
    return a / (a + b);
}

inline void die(const std::string &epitaph, int code=1)
{
    fputs((epitaph + '\n').c_str(), stderr);
    soft_exit(code);
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
// You'd think this would be in the C library.
// Like strchr, but with a max length so it doesn't
// run over the end of the buffer.  Basically,
// strings suck in C.
//

    inline char *
strnchr(char *str, char charToFind, size_t maxLen)
{
    for (size_t i = 0; i < maxLen; i++) {
        if (str[i] == charToFind) {
            return str + i;
        }
        if (str[i] == 0) {
            return NULL;
        }
    }
    return NULL;
}

    inline const char *
strnchr(const char *str, char charToFind, size_t maxLen)
{
    for (size_t i = 0; i < maxLen; i++) {
        if (str[i] == charToFind) {
            return str + i;
        }
        if (str[i] == 0) {
            return NULL;
        }
    }
    return NULL;
}
    
// Check whether a string str ends with a given pattern
    inline bool
stringEndsWith(const char* str, const char* pattern)
{
    if (strlen(str) < strlen(pattern)) {
        return false;
    } else {
#ifdef _MSC_VER
        return _stricmp(str + (strlen(str) - strlen(pattern)), pattern) == 0;
#else
        return strcmp(str + (strlen(str) - strlen(pattern)), pattern) == 0;
#endif
    }
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
    
    // from MurmurHash3, public domain, http://code.google.com/p/smhasher/wiki/MurmurHash3

#ifdef _MSC_VER
#define ROTL32(x,y)     _rotl(x,y)
#define ROTL64(x,y)     _rotl64(x,y)
#define BIG_CONSTANT(x) (x)
#else
#define ROTL32(x,y)     rotl32(x,y)
#define ROTL64(x,y)     rotl64(x,y)
#define BIG_CONSTANT(x) (x##LLU)
#endif

    inline _uint64
fmix64(
    _uint64 k)
 
{
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;
  return k;
}
    
    inline _uint64
hash64(
    const void* key,
    int len)
{
    const _uint8* data = (const _uint8*) key;
    const int nblocks = len / 16;

    _uint64 h1 = 0x460a5856818aaba3LL;
    _uint64 h2 = h1;

    const _uint64 c1 = BIG_CONSTANT(0x87c37b91114253d5);
    const _uint64 c2 = BIG_CONSTANT(0x4cf5ad432745937f);

    //----------
    // body

    const _uint64 * blocks = (const _uint64 *)(data);

    for(int i = 0; i < nblocks; i++)
    {
      _uint64 k1 = blocks[i*2+0];
      _uint64 k2 = blocks[i*2+1];

      k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;

      h1 = ROTL64(h1,27); h1 += h2; h1 = h1*5+0x52dce729;

      k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

      h2 = ROTL64(h2,31); h2 += h1; h2 = h2*5+0x38495ab5;
    }

    //----------
    // tail

    const _uint8 * tail = (const _uint8*)(data + nblocks*16);

    _uint64 k1 = 0;
    _uint64 k2 = 0;

    switch(len & 15)
    {
    case 15: k2 ^= ((_uint64)tail[14]) << 48;
    case 14: k2 ^= ((_uint64)tail[13]) << 40;
    case 13: k2 ^= ((_uint64)tail[12]) << 32;
    case 12: k2 ^= ((_uint64)tail[11]) << 24;
    case 11: k2 ^= ((_uint64)tail[10]) << 16;
    case 10: k2 ^= ((_uint64)tail[ 9]) << 8;
    case  9: k2 ^= ((_uint64)tail[ 8]) << 0;
             k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

    case  8: k1 ^= ((_uint64)tail[ 7]) << 56;
    case  7: k1 ^= ((_uint64)tail[ 6]) << 48;
    case  6: k1 ^= ((_uint64)tail[ 5]) << 40;
    case  5: k1 ^= ((_uint64)tail[ 4]) << 32;
    case  4: k1 ^= ((_uint64)tail[ 3]) << 24;
    case  3: k1 ^= ((_uint64)tail[ 2]) << 16;
    case  2: k1 ^= ((_uint64)tail[ 1]) << 8;
    case  1: k1 ^= ((_uint64)tail[ 0]) << 0;
             k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
    };

    //----------
    // finalization

    h1 ^= len; h2 ^= len;

    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    return h2;
}

} // namespace util

int FirstPowerOf2GreaterThanOrEqualTo(int value);

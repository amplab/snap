#pragma once

#include <map>
#include "stdafx.h"
#include "Compat.h"
#include "Tables.h"
#include "exit.h"
#include "GenericFile.h"
using std::max;
using std::min;

//
// General utilities.
//

namespace util {

#define BEGEND(container) (container).begin(), (container).end()

inline double ratio(double a, double b=1)
{
    return a / (a + b);
}

//
// Turn the value into a string with comma formatting (so 1,234,567 instead of 1234567).
// Produces a null-terminated string.
//
extern char *FormatUIntWithCommas(_uint64 val, char *outputBuffer, size_t outputBufferSize);

const int MAXLINE = 1024;

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
    if (bases != NULL && bases != rc) {
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

    inline void
diffSequence(
    char* io_sequence,
    int length,
    const char* reference,
    int refLength,
    int offset) // offset of sequence from reference
{
    int overlap = min(offset + length, refLength) - max(offset, 0); // # bases offset
    int start = max(-offset, 0); // index in sequence
    int last = 0;
    for (int i = 0; i < overlap; i++) {
        int is = start + i;
        int ir = start + i + offset;
        int found = 0;
        for (int j = 0; j < 3; j++) {
            int delta = j == 0 ? last :
                j == 1 ? (last ? 0 : -1) :
                (last ? -last : +1);
            if (ir + delta >= 0 && ir + delta < refLength && toupper(reference[ir + delta]) == toupper(io_sequence[is])) {
                io_sequence[is] = (islower(io_sequence[is]) ? "-,+" : "-.+")[1 + delta];
                found = delta;
                break;
            }
        }
        last = found;
    }
}

    inline void
toLower(
    char* buffer,
    int length)
{
    for (int i = 0; i < length; i++) {
        buffer[i] = tolower(buffer[i]);
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
    
    inline int
log10bucket(
    int n)
{
    unsigned factor = 1;
    if (n < 0) {
        n = -n;
        factor = -1;
    }
    while (n >= 10) {
        n /= 10;
        factor *= 10;
    }
    return n * factor;
}
    
    inline _int64
log10bucket(
    _int64 n)
{
    _int64 factor = 1;
    while (abs(n) >= 10) {
        n /= 10;
        factor *= 10;
    }
    return n * factor;
}
    // from MurmurHash3, public domain, http://code.google.com/p/smhasher/wiki/MurmurHash3
#define ROTL32(x,y)     bit_rotate_left(x,y)
#define ROTL64(x,y)     bit_rotate_left64(x,y)

#ifdef _MSC_VER
#define BIG_CONSTANT(x) (x)
#else
#define BIG_CONSTANT(x) (x##LLU)
#endif

    inline _uint32
fmix32(
    _uint32 h)
{
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

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
    inline _uint32
hash(
    const void* key,
    int len)
{
    const _uint8* data = (const _uint8*) key;
    const int nblocks = len / 4;
 
    _uint32 h1 = 0x811f6d67; // seed, const from a random guid for now
 
    const _uint32 c1 = 0xcc9e2d51;
    const _uint32 c2 = 0x1b873593;
 
    //----------
    // body
 
    const _uint32 * blocks = (const _uint32 *)(data + nblocks*4);
 
    for(int i = -nblocks; i; i++)
    {
        _uint32 k1 = blocks[i];
 
        k1 *= c1;
        k1 = ROTL32(k1,15);
        k1 *= c2;
    
        h1 ^= k1;
        h1 = ROTL32(h1,13); 
        h1 = h1*5+0xe6546b64;
    }
 
    //----------
    // tail
 
    const _uint32 * tail = (const _uint32*)(data + nblocks*4);
    _uint32 k1 = 0;
 
    switch(len & 3)
    {
        case 3: k1 ^= tail[2] << 16;
        case 2: k1 ^= tail[1] << 8;
        case 1: k1 ^= tail[0];
                k1 *= c1; k1 = ROTL32(k1,15); k1 *= c2; h1 ^= k1;
    };
 
    //----------
    // finalization
 
    h1 ^= len;
    h1 = fmix32(h1);
    return h1;
}
    
    inline _uint64
hash64(
    _uint64 x)
{
    return fmix64(x);
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

struct IdPair
{
    unsigned id, value;
    bool operator==(const IdPair& b) const
    {
        return id == b.id && value == b.value;
    }
    static bool comparator(const IdPair& a, const IdPair& b)
    {
        return a.id < b.id || (a.id == b.id && a.value < b.value);
    }
    static bool valueComparator(const IdPair& a, const IdPair& b)
    {
        return a.value < b.value || (a.value == b.value && a.id < b.id);
    }
    static bool valueComparatorDescending(const IdPair& a, const IdPair& b)
    {
        return a.value > b.value;
    }
    IdPair() : id(0), value(0) {}
    IdPair(unsigned i_id, unsigned i_value) : id(i_id), value(i_value){}
    // for use as key in VariableSizeMap
    IdPair(int i) : id((unsigned) i), value(0) {}
    bool operator==(int x) const
    { return id == (unsigned) x && value == 0; }
    bool operator!=(int x) const
    { return id != (unsigned) x || value != 0; }
    operator _uint64()
    { return (((_uint64) id) << 32) | (_uint32) value; }
};
    
struct IdIntPair
{
    unsigned id;
    int value;
    bool operator==(const IdIntPair& b) const
    {
        return id == b.id && value == b.value;
    }
    static bool comparator(const IdIntPair& a, const IdIntPair& b)
    {
        return a.id < b.id || (a.id == b.id && a.value < b.value);
    }
    static bool valueComparator(const IdIntPair& a, const IdIntPair& b)
    {
        return a.value < b.value;
    }
    IdIntPair()
        : id(0), value(0)
    {}
    IdIntPair(unsigned i_id, int i_value)
        : id(i_id), value(i_value)
    {}
    IdIntPair(_uint64 x)
        : id((unsigned) (x >> 32)), value((int) x)
    {}
    // for use as key in VariableSizeMap
    IdIntPair(int i) : id((unsigned) i), value(0) {}
    bool operator==(int x) const
    { return id == (unsigned) x && value == 0; }
    bool operator!=(int x) const
    { return id != (unsigned) x || value != 0; }
    operator _uint64()
    { return (((_uint64) id) << 32) | (_uint32) value; }
};


void memrevcpy(void* dst, const void* src, size_t bytes);


} // namespace util

_int64 FirstPowerOf2GreaterThanOrEqualTo(_int64 value);
int cheezyLogBase2(_int64 value);

//
// Check if a is within distance of b, coping properly with the varagies of unsigneds.
// There's a similar function for GenomeLocations defined in Genome.h.
//
inline bool isWithin(unsigned a, unsigned b, unsigned distance)
    {
	return a <= b && a+distance >= b || a >= b && a <= b + distance;
}

inline int getSignBit64(_int64 value)
{
    return (value >> 63) & 1;
}

inline int getSignBit32(int value)
{
    return (value >> 31) & 1;
}


// Utility class for synchronization: NWaiter.
// This class is initialized with a number, n.
// It has two public methods: wait() and signal().
// wait() will block until signal() has been called n times.
class NWaiter
{
public:
	void wait();
	void signal();
	NWaiter(size_t n);
	~NWaiter();

private:
	size_t _signalsRequired;
	size_t _signalsReceived;
	EventObject _waiter;
	ExclusiveLock _lock;
};

//
// Version of fgets that dynamically (re-)allocates the buffer to be big enough to fit the whole line
//
char *reallocatingFgets(char **buffer, int *io_bufferSize, FILE *stream);
char *reallocatingFgetsGenericFile(char **buffer, int *io_bufferSize, GenericFile *file);



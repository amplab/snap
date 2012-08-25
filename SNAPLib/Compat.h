/*++

Module Name:

    compat.h

Abstract:

    Functions that provide compatibility between the Windows and Linux versions,
    and mostly that serve to keep #ifdef's out of the main code in order to
    improve readibility.

Authors:

    Bill Bolosky, November, 2011

Environment:

    User mode service.

Revision History:


--*/

#pragma once

#ifdef  _MSC_VER

typedef unsigned _int64 _uint64;
typedef unsigned _int32 _uint32;

// <http://stackoverflow.com/questions/126279/c99-stdint-h-header-and-ms-visual-studio>
const _int64 INT64_MAX = MAXINT64;
const _int64 INT64_MIN = MININT64;
const _uint32 UINT32_MAX = MAXUINT32;
const _int32 INT32_MIN = MININT32;
const _int32 INT32_MAX = MAXINT32;

static const double LOG10 = log(10.0);
inline double exp10(double x) { return exp(x * LOG10); }

const void* memmem(const void* data, const size_t dataLength, const void* pattern, const size_t patternLength);

typedef CRITICAL_SECTION    ExclusiveLock;
typedef HANDLE SingleWaiterObject;      // This is an event in Windows.  It's just a synchronization object that you can wait for and set.

#define PATH_SEP '\\'
#define snprintf _snprintf
#define mkdir(path, mode) _mkdir(path)

// <http://stackoverflow.com/questions/9021502/whats-the-difference-between-strtok-r-and-strtok-s-in-c>
#define strtok_r strtok_s
#define strncasecmp _strnicmp
#define atoll(S) _atoi64(S)


#else   // _MSC_VER

#include <pthread.h>
// <http://stackoverflow.com/questions/986426/what-do-stdc-limit-macros-and-stdc-constant-macros-mean>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <assert.h>
#include <float.h>

#ifdef __linux__
#include <sched.h>  // For sched_setaffinity
#endif

typedef int64_t _int64;
typedef uint64_t _uint64;
typedef int32_t _int32;
typedef uint32_t _uint32;
typedef uint8_t BYTE;
typedef void *PVOID;

// TODO: check if Linux libs have exp10 function
#include <math.h>
static const double LOG10 = log(10.0);
inline double exp10(double x) { return exp(x * LOG10); }

#define __in /* nothing */

#define PATH_SEP '/'

#ifdef DEBUG
#define _ASSERT assert
#else
#define _ASSERT(x) {}
#endif

#define __min(x,y) ((x)<(y) ? (x) : (y))
#define __max(x,y) ((x)>(y) ? (x) : (y))
#define MAX_PATH 4096
#define __cdecl __attribute__((__cdecl__))

inline bool _BitScanForward64(unsigned long *result, _uint64 x) {
    *result = __builtin_ctzll(x);
    return x != 0;
}

// We implement SingleWaiterObject using a mutex because POSIX unnamed semaphores don't work on OS X
class SingleWaiterObjectImpl;

typedef pthread_mutex_t ExclusiveLock;
typedef SingleWaiterObjectImpl *SingleWaiterObject;

#endif  // _MSC_VER

//
// Get the time since some predefined time.  The predefined time must not change during any particular program run.
//
_int64 timeInMillis();
_int64 timeInNanos();

//
// Exclusive locks.  These have the obvious semantics: At most one thread can acquire one at any time, the others block
// until the first one releases it.
//
void AcquireExclusiveLock(ExclusiveLock *lock);
void ReleaseExclusiveLock(ExclusiveLock *lock);
bool InitializeExclusiveLock(ExclusiveLock *lock);
bool DestroyExclusiveLock(ExclusiveLock *lock);

//
// Single waiter objects.  The semantics are that a single thread can wait on one of these, and when it's
// set by any thread, the waiter will proceed.  It works regardless of the order of waiting and signalling.
// Can be reset back to unsignalled state.
//


bool CreateSingleWaiterObject(SingleWaiterObject *newWaiter);
void DestroySingleWaiterObject(SingleWaiterObject *waiter);
void SignalSingleWaiterObject(SingleWaiterObject *singleWaiterObject);
bool WaitForSingleWaiterObject(SingleWaiterObject *singleWaiterObject);
void ResetSingleWaiterObject(SingleWaiterObject *singleWaiterObject);


//
// Thread-safe read-modify-write operations
//
_uint32 InterlockedIncrementAndReturnNewValue(volatile _uint32 *valueToIncrement);
int InterlockedDecrementAndReturnNewValue(volatile int *valueToDecrement);
_int64 InterlockedAdd64AndReturnNewValue(volatile _int64 *valueToWhichToAdd, _int64 amountToAdd);
_uint32 InterlockedCompareExchange32AndReturnOldValue(volatile _uint32 *valueToUpdate, _uint32 replacementValue, _uint32 desiredPreviousValue);
_uint64 InterlockedCompareExchange64AndReturnOldValue(volatile _uint64 *valueToUpdate, _uint64 replacementValue, _uint64 desiredPreviousValue);
void* InterlockedCompareExchangePointerAndReturnOldValue(volatile void **valueToUpdate, void* replacementValue, void* desiredPreviousValue);

//
// Functions for creating and binding threads.
//
typedef void (*ThreadMainFunction) (void *threadMainFunctionParameter);
bool StartNewThread(ThreadMainFunction threadMainFunction, void *threadMainFunctionParameter);
void BindThreadToProcessor(unsigned processorNumber); // This hard binds a thread to a processor.  You can no-op it at some perf hit.

unsigned GetNumberOfProcessors();

_int64 QueryFileSize(const char *fileName);

// returns true on success
bool DeleteSingleFile(const char* filename); // DeleteFile is a Windows macro...

class LargeFileHandle;

// open binary file, supports "r" for read, "w" for rewrite/create, "a" for append
LargeFileHandle* OpenLargeFile(const char* filename, const char* mode);

size_t WriteLargeFile(LargeFileHandle* file, void* buffer, size_t bytes);

size_t ReadLargeFile(LargeFileHandle* file, void* buffer, size_t bytes);

// closes and deallocates
void CloseLargeFile(LargeFileHandle* file);

// open and close memory mapped files
// currently just readonly, could add flags for r/w if necessary

class MemoryMappedFile;

MemoryMappedFile* OpenMemoryMappedFile(const char* filename, size_t offset, size_t length, void** o_contents, bool write = false, bool sequential = false);

// closes and deallocates the file structure
void CloseMemoryMappedFile(MemoryMappedFile* mappedFile);

class AsyncFile
{
public:

    // open a new file for reading and/or writing
    static AsyncFile* open(const char* filename, bool write);

    // free resources; must have destroyed all readers & writers first
    virtual bool close() = 0;

    // abstract class for asynchronous writes
    class Writer
    {
    public:
        // waits for all writes to complete, frees resources
        virtual bool close() = 0;

        // begin a write; if there is already a write in progress, might wait for it to complete
        virtual bool beginWrite(void* buffer, size_t length, size_t offset, size_t *bytesWritten) = 0;

        // wait for all prior beginWrites to complete
        virtual bool waitForCompletion() = 0;
    };

    // get a new writer, e.g. for another thread to use
    virtual Writer* getWriter() = 0;

    // abstract class for asynchronous reads
    class Reader
    {
    public:
        // waits for alls reads to complete, frees resources
        virtual bool close() = 0;

        // begin a new read; if there is already a read in progress, might wait for it to complete
        virtual bool beginRead(void* buffer, size_t length, size_t offset, size_t *bytesRead) = 0;

        // wait for all prior beginReads to complete
        virtual bool waitForCompletion() = 0;
    };

    // get a new reader, e.g. for another thread to use
    virtual Reader* getReader() = 0;
};

//
// Macro for counting trailing zeros of a 64-bit value
//
#ifdef _MSC_VER
#define CountTrailingZeroes(x, ans) {_BitScanForward64(&ans, x);}
#else
#define CountTrailingZeroes(x, ans) {ans = __builtin_ctzll(x);}
#endif

//
// 64 bit version of fseek.
//

int _fseek64bit(FILE *stream, _int64 offset, int origin);

#ifndef _MSC_VER

#define MININT32   ((int32_t)  0x80000000)
#define MAXINT32   ((int32_t)  0x7fffffff)

#endif

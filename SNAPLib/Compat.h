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
#include <Windows.h>

typedef unsigned _int64 _uint64;
typedef unsigned _int32 _uint32;
typedef unsigned char _uint8;
typedef unsigned short _uint16;

// <http://stackoverflow.com/questions/126279/c99-stdint-h-header-and-ms-visual-studio>
const _uint64 UINT64_MAX = MAXUINT64;
const _int64 INT64_MAX = MAXINT64;
const _int64 INT64_MIN = MININT64;
const _uint32 UINT32_MAX = MAXUINT32;
const _int32 INT32_MIN = MININT32;
const _int32 INT32_MAX = MAXINT32;
const _uint16 UINT16_MAX = MAXUINT16;
const _int16 INT16_MAX = MAXINT16;
const _int16 INT16_MIN = MININT16;

static const double LOG10 = log(10.0);
inline double exp10(double x) { return exp(x * LOG10); }

const void* memmem(const void* data, const size_t dataLength, const void* pattern, const size_t patternLength);

typedef CRITICAL_SECTION    UnderlyingExclusiveLock;
typedef HANDLE SingleWaiterObject;      // This is an event in Windows.  It's just a synchronization object that you can wait for and set.  "Single" means only one thread can wait on it at a time.
typedef HANDLE EventObject;

#define PATH_SEP '\\'
#define snprintf _snprintf
#define mkdir(path, mode) _mkdir(path)
#define strdup(s) _strdup(s)

// <http://stackoverflow.com/questions/9021502/whats-the-difference-between-strtok-r-and-strtok-s-in-c>
#define strtok_r strtok_s
#define strncasecmp _strnicmp
#define atoll(S) _atoi64(S)

#define bit_rotate_right(value, shift) _rotr(value, shift)
#define bit_rotate_left(value, shift) _rotl(value, shift)
#define bit_rotate_right64(value, shift) _rotr64(value, shift)
#define bit_rotate_left64(value, shift) _rotl64(value, shift)

int getpagesize();
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

#ifndef __APPLE__
#include <xmmintrin.h>  // This is currently (in Dec 2013) broken on Mac OS X 10.9 (Apple clang-500.2.79)
#else
#define _mm_prefetch(...) {}
#endif

typedef int64_t _int64;
typedef uint64_t _uint64;
typedef int32_t _int32;
typedef uint32_t _uint32;
typedef uint16_t _uint16;
typedef int16_t _int16;
typedef uint8_t BYTE;
typedef uint8_t _uint8;
typedef int8_t _int8;
typedef void *PVOID;

// TODO: check if Linux libs have exp10 function
#include <math.h>
static const double LOG10 = log(10.0);
inline double exp10(double x) { return exp(x * LOG10); }

#define __in /* nothing */

#define PATH_SEP '/'

#ifdef DEBUG
#define _ASSERT assert
#ifndef _DEBUG
#define _DEBUG 1	// Compat with Windows version
#endif // !_DEBUG
#else
#define _ASSERT(x) {}
#endif

#define __min(x,y) ((x)<(y) ? (x) : (y))
#define __max(x,y) ((x)>(y) ? (x) : (y))
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#define MAX_PATH 4096
#define __cdecl __attribute__((__cdecl__))

#define _stricmp strcasecmp

inline bool _BitScanForward64(unsigned long *result, _uint64 x) {
    *result = __builtin_ctzll(x);
    return x != 0;
}

// We implement SingleWaiterObject using a mutex because POSIX unnamed semaphores don't work on OS X
class SingleWaiterObjectImpl;

typedef pthread_mutex_t UnderlyingExclusiveLock;
typedef SingleWaiterObjectImpl *SingleWaiterObject; // "Single" means only one thread can wait on it at a time.

class EventObjectImpl;
typedef EventObjectImpl *EventObject;

inline unsigned bit_rotate_right(unsigned value, unsigned shift)
{
    if (shift%32 == 0) return value;
    return value >> (shift%32) | (value << (32 - shift%32));
}

inline unsigned bit_rotate_left(unsigned value, unsigned shift)
{
    if (shift%32 == 0) return value;
    return value << (shift %32) | (value >> (32 - shift%32));
}

inline _uint64 bit_rotate_right64(_uint64 value, unsigned shift)
{
    if (shift%64 == 0) return value;
    return value >> (shift%64) | (value << (64 - shift%64));
}

inline _uint64 bit_rotate_left64(_uint64 value, unsigned shift)
{
    if (shift%64 == 0) return value;
    return value << (shift %64) | (value >> (64 - shift%64));
}

#endif  // _MSC_VER

struct NamedPipe;	// It's bi-directional, which in Unix means it's actually two pipes

extern NamedPipe *OpenNamedPipe(const char *pipeName, bool serverSide);
extern bool ReadFromNamedPipe(NamedPipe *pipe, char *outputBuffer, size_t outputBufferSize);
extern bool WriteToNamedPipe(NamedPipe *pipe, const char *stringToWrite);	// Null-terminated string
extern void CloseNamedPipe(NamedPipe *pipe);

extern const char *DEFAULT_NAMED_PIPE_NAME;

//
// Get the time since some predefined time.  The predefined time must not change during any particular program run.
//
_int64 timeInMillis();
_int64 timeInNanos();

//#define PROFILE_WAIT

void PrintWaitProfile();


//
// Exclusive locks.  These have the obvious semantics: At most one thread can acquire one at any time, the others block
// until the first one releases it.  In the DEBUG build we wrap the lock in a class that ensures that it's initialized before
// it's used (which we found out the hard way isn't always so obvious).
//
extern void AcquireUnderlyingExclusiveLock(UnderlyingExclusiveLock *);
bool InitializeUnderlyingExclusiveLock(UnderlyingExclusiveLock *lock);
void ReleaseUnderlyingExclusiveLock(UnderlyingExclusiveLock *lock);
bool DestroyUnderlyingExclusiveLock(UnderlyingExclusiveLock *lock);

#ifdef _DEBUG
class ExclusiveLock {
public:
    UnderlyingExclusiveLock lock;
    bool                    initialized;
	bool					wholeProgramScope;

#ifdef _MSC_VER
    DWORD                   holderThreadId;
#endif // _MSC_VER


    ExclusiveLock() : initialized(false), holderThreadId(0), wholeProgramScope(false) {}
    ~ExclusiveLock() {_ASSERT(!initialized || wholeProgramScope);}   // Must DestroyExclusiveLock first
};

inline void SetExclusiveLockWholeProgramScope(ExclusiveLock *lock)
{
	lock->wholeProgramScope = true;
}

inline void AcquireExclusiveLock(ExclusiveLock *lock)
{
    _ASSERT(lock->initialized);
    AcquireUnderlyingExclusiveLock(&lock->lock);
#ifdef _MSC_VER
    // If you see this go off, you're probably trying a recursive lock acquisition (i.e., twice on the same thead), 
    // which is legal in Windows and a deadlock in Linux.
    _ASSERT(lock->holderThreadId == 0);                
    lock->holderThreadId = GetCurrentThreadId();
#endif // _MSC_VER

}

inline void AssertExclusiveLockHeld(ExclusiveLock *lock)
{
#ifdef _MSC_VER
    _ASSERT(GetCurrentThreadId() == lock->holderThreadId);
#endif // _MSC_VER

}

inline bool InitializeExclusiveLock(ExclusiveLock *lock)
{
    _ASSERT(!lock->initialized);
    lock->initialized = true;
    return InitializeUnderlyingExclusiveLock(&lock->lock);
}

inline void ReleaseExclusiveLock (ExclusiveLock *lock)
{
    _ASSERT(lock->initialized);
#ifdef _MSC_VER
    _ASSERT(GetCurrentThreadId() == lock->holderThreadId);
    lock->holderThreadId = 0;
#endif // _MSC_VER

    ReleaseUnderlyingExclusiveLock(&lock->lock);
}

inline void DestroyExclusiveLock(ExclusiveLock *lock)
{
#ifdef _MSC_VER
    _ASSERT(lock->holderThreadId == 0);
#endif // _MSC_VER

	_ASSERT(!lock->wholeProgramScope);
    _ASSERT(lock->initialized);
    lock->initialized = false;
    DestroyUnderlyingExclusiveLock(&lock->lock);
}
#else   // _DEBUG
#define ExclusiveLock UnderlyingExclusiveLock
#define InitializeExclusiveLock InitializeUnderlyingExclusiveLock
#define ReleaseExclusiveLock ReleaseUnderlyingExclusiveLock
#define DestroyExclusiveLock DestroyUnderlyingExclusiveLock
#define AssertExclusiveLockHeld(l) /* nothing */
#define SetExclusiveLockWholeProgramScope(l) /*nothing*/
#endif // _DEBUG

#ifdef PROFILE_WAIT
#define AcquireExclusiveLock(lock) AcquireExclusiveLockProfile((lock), __FUNCTION__, __LINE__)
void AcquireExclusiveLockProfile(ExclusiveLock *lock, const char* fn, int line);
#elif   _DEBUG
// already defined above
#else   // !debug, !profile_wait
#define AcquireExclusiveLock AcquireUnderlyingExclusiveLock
#endif


//
// Single waiter objects.  The semantics are that a single thread can wait on one of these, and when it's
// set by any thread, the waiter will proceed.  It works regardless of the order of waiting and signalling.
// Can be reset back to unsignalled state.
//


bool CreateSingleWaiterObject(SingleWaiterObject *newWaiter);
void DestroySingleWaiterObject(SingleWaiterObject *waiter);
void SignalSingleWaiterObject(SingleWaiterObject *singleWaiterObject);
#ifdef PROFILE_WAIT
#define WaitForSingleWaiterObject(o) WaitForSingleWaiterObjectProfile((o), __FUNCTION__, __LINE__)
bool WaitForSingleWaiterObjectProfile(SingleWaiterObject *singleWaiterObject, const char* fn, int line);
#else
bool WaitForSingleWaiterObject(SingleWaiterObject *singleWaiterObject);
#endif
void ResetSingleWaiterObject(SingleWaiterObject *singleWaiterObject);

//
// An Event is a synchronization object that acts as a gateway: it can either be open
// or closed.  Open events allow all waiters to proceed, while closed ones block all
// waiters.  Events can be opened and closed multiple times, and can have any number of
// waiters.
//

void CreateEventObject(EventObject *newEvent);
void DestroyEventObject(EventObject *eventObject);
void AllowEventWaitersToProceed(EventObject *eventObject);
void PreventEventWaitersFromProceeding(EventObject *eventObject);
#ifdef PROFILE_WAIT
#define WaitForEvent(o) WaitForEventProfile((o), __FUNCTION__, __LINE__)
void WaitForEventProfile(EventObject *eventObject, const char* fn, int line);
#else
void WaitForEvent(EventObject *eventObject);
#endif
bool WaitForEventWithTimeout(EventObject *eventObject, _int64 timeoutInMillis); // Returns true if the event was set, false if the timeout happened

//
// Thread-safe read-modify-write operations
//
int InterlockedIncrementAndReturnNewValue(volatile int *valueToIncrement);
int InterlockedDecrementAndReturnNewValue(volatile int *valueToDecrement);
_int64 InterlockedAdd64AndReturnNewValue(volatile _int64 *valueToWhichToAdd, _int64 amountToAdd);
_uint32 InterlockedCompareExchange32AndReturnOldValue(volatile _uint32 *valueToUpdate, _uint32 replacementValue, _uint32 desiredPreviousValue);
_uint64 InterlockedCompareExchange64AndReturnOldValue(volatile _uint64 *valueToUpdate, _uint64 replacementValue, _uint64 desiredPreviousValue);
void* InterlockedCompareExchangePointerAndReturnOldValue(void * volatile *valueToUpdate, void* replacementValue, void* desiredPreviousValue);

//
// Functions for creating and binding threads.
//
typedef void (*ThreadMainFunction) (void *threadMainFunctionParameter);
bool StartNewThread(ThreadMainFunction threadMainFunction, void *threadMainFunctionParameter);
void BindThreadToProcessor(unsigned processorNumber); // This hard binds a thread to a processor.  You can no-op it at some perf hit.
#ifdef  _MSC_VER
#define GetThreadId() GetCurrentThreadId()
#else   // _MSC_VER
#define GetThreadId() pthread_self()
#endif  // _MSC_VER

void SleepForMillis(unsigned millis);

unsigned GetNumberOfProcessors();

_int64 QueryFileSize(const char *fileName);

// returns true on success
bool DeleteSingleFile(const char* filename); // DeleteFile is a Windows macro...

// returns true on success
bool MoveSingleFile(const char* oldFileName, const char* newFileName);

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

//
void AdviseMemoryMappedFilePrefetch(const MemoryMappedFile *mappedFile);

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
#define CountLeadingZeroes(x, ans) {_BitScanReverse64(&ans, x);}
#define CountTrailingZeroes(x, ans) {_BitScanForward64(&ans, x);}
#define ByteSwapUI64(x) (_byteswap_uint64(x))
#else
#define CountLeadingZeroes(x, ans) {ans = __builtin_clzll(x);}
#define CountTrailingZeroes(x, ans) {ans = __builtin_ctzll(x);}
#define ByteSwapUI64(x) (__builtin_bswap64(x))
#endif

//
// 64 bit version of fseek.
//

int _fseek64bit(FILE *stream, _int64 offset, int origin);

#ifndef _MSC_VER

#define MININT32   ((int32_t)  0x80000000)
#define MAXINT32   ((int32_t)  0x7fffffff)

#endif

//
// Class for handling mapped files.  It's got the same interface for both platforms, but different implementations.
//
class FileMapper {
public:
    FileMapper();
    ~FileMapper();

    // can only be called once - only usable for a single file
    bool init(const char *fileName);

    const size_t getFileSize() {
        _ASSERT(initialized);
        return fileSize;
    }

    // can get multiple mappings on the same file
    char *createMapping(size_t offset, size_t amountToMap, void** o_token);

    // MUST call unmap on each token out of createMapping, the destructor WILL NOT cleanup
    void unmap(void* token);

private:
    bool        initialized;
    const char* fileName;
    size_t      fileSize;
    size_t      pagesize;
    int         mapCount; // simple count of mappings that have not yet been unmapped

#ifdef  _MSC_VER
    HANDLE      hFile;
    HANDLE      hMapping;

    _int64 millisSpentInReadFile;
    _int64 countOfImmediateCompletions;
    _int64 countOfDelayedCompletions;
    _int64 countOfFailures;
#else   // _MSC_VER
    static const int madviseSize = 4 * 1024 * 1024;
    typedef std::pair<void*,size_t> UnmapToken;

    int         fd;
    _uint64     lastPosMadvised;
#endif  // _MSC_VER
};

//
// Call to keep the OS from putting the machine asleep
//
void PreventMachineHibernationWhileThisThreadIsAlive();

//
// Reduce our scheduling priority to be nicer to other jobs.
//
void SetToLowSchedulingPriority();


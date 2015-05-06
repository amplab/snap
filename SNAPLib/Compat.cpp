/*++

Module Name:

    compat.cpp

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
#include "stdafx.h"
#include "Compat.h"
#include "BigAlloc.h"
#ifndef _MSC_VER
#include <fcntl.h>
#include <aio.h>
#include <err.h>
#include <unistd.h>
#include <signal.h>
#endif
#include "exit.h"
#ifdef PROFILE_WAIT
#include <map>
#endif
#include "DataWriter.h"
#include "Error.h"

using std::min;
using std::max;

#ifdef PROFILE_WAIT

#undef AcquireExclusiveLock
#undef WaitForSingleWaiterObject
#undef WaitForEvent

void AcquireUnderlyingExclusiveLock(UnderlyingExclusiveLock *lock);
bool WaitForSingleWaiterObject(SingleWaiterObject *singleWaiterObject);
void WaitForEvent(EventObject *eventObject); 

using std::map;
using std::string;
static map<string,_int64> times;

void addTime(const char* fn, int line, _int64 time)
{
    if (time > 0) {
        char s[20];
        sprintf(s, ":%d", line);
        string key = string(fn) + string(s);
        times[key] += time;
    }
}

void AcquireExclusiveLockProfile(ExclusiveLock *lock, const char* fn, int line)
{
    _int64 start = timeInMillis();
    AcquireExclusiveLock(lock);
    addTime(fn, line, timeInMillis() - start);
}

bool WaitForSingleWaiterObjectProfile(SingleWaiterObject *singleWaiterObject, const char* fn, int line)
{
    _int64 start = timeInMillis();
    bool result = WaitForSingleWaiterObject(singleWaiterObject);
    addTime(fn, line, timeInMillis() - start);
    return result;
}

void WaitForEventProfile(EventObject *eventObject, const char* fn, int line)
{
    _int64 start = timeInMillis();
    WaitForEvent(eventObject); 
    addTime(fn, line, timeInMillis() - start);
}

#endif

void PrintWaitProfile()
{
#ifdef PROFILE_WAIT
    printf("function:line    wait_time (s)\n");
    for (map<string,_int64>::iterator lt = times.begin(); lt != times.end(); lt++) {
        printf("%s %.3f\n", lt->first.data(), lt->second * 0.0001);
    }
#endif
}

#ifdef  _MSC_VER

const void* memmem(const void* data, const size_t dataLength, const void* pattern, const size_t patternLength)
{
    if (dataLength < patternLength) {
        return NULL;
    }
    const void* p = data;
    const char first = *(char*)pattern;
    size_t count = dataLength - patternLength + 1;
    while (count > 0) {
        const void* q = memchr(p, first, count);
        if (q == NULL) {
            return NULL;
        }
        if (0 == memcmp(q, pattern, patternLength)) {
            return q;
        }
        count -= ((char*)q - (char*)p) + 1;
        p = (char*)q + 1;
    }
    return NULL;
}

_int64 timeInMillis()
/**
 * Get the current time in milliseconds since some arbitrary starting point
 * (e.g. system boot or epoch)
 */
{
    return GetTickCount64();
}

_int64 timeInNanos()
{
    static _int64 performanceFrequency = 0;

    if (0 == performanceFrequency) {
        QueryPerformanceFrequency((PLARGE_INTEGER)&performanceFrequency);
    }

    LARGE_INTEGER perfCount;
    QueryPerformanceCounter(&perfCount);
    return perfCount.QuadPart * 1000000000 / performanceFrequency;
}

void AcquireUnderlyingExclusiveLock(UnderlyingExclusiveLock *lock) {
    EnterCriticalSection(lock);
}

void ReleaseUnderlyingExclusiveLock(UnderlyingExclusiveLock *lock) {
    LeaveCriticalSection(lock);
}

bool InitializeUnderlyingExclusiveLock(UnderlyingExclusiveLock *lock) {
    InitializeCriticalSection(lock);
    return true;
}

bool DestroyUnderlyingExclusiveLock(UnderlyingExclusiveLock *lock) {
    DeleteCriticalSection(lock);
    return true;
}

bool CreateSingleWaiterObject(SingleWaiterObject *newWaiter)
{
    *newWaiter = CreateEvent(NULL,TRUE,FALSE,NULL);
    if (NULL == *newWaiter) {
        return false;
    }
    return true;
}

void DestroySingleWaiterObject(SingleWaiterObject *waiter)
{
    CloseHandle(*waiter);
}

void SignalSingleWaiterObject(SingleWaiterObject *singleWaiterObject) {
    SetEvent(*singleWaiterObject);
}

bool WaitForSingleWaiterObject(SingleWaiterObject *singleWaiterObject) {
    DWORD retVal = WaitForSingleObject(*singleWaiterObject,INFINITE);
    if (WAIT_OBJECT_0 != retVal) {
        return false;
    }
    return true;
}

void ResetSingleWaiterObject(SingleWaiterObject *singleWaiterObject) {
    ResetEvent(*singleWaiterObject);
}

//
// In Windows, the single waiter objects are already implemented using events.
//
void CreateEventObject(EventObject *newEvent) {CreateSingleWaiterObject(newEvent);}
void DestroyEventObject(EventObject *eventObject) {DestroySingleWaiterObject(eventObject);}
void AllowEventWaitersToProceed(EventObject *eventObject) {SignalSingleWaiterObject(eventObject);}
void PreventEventWaitersFromProceeding(EventObject *eventObject) {ResetSingleWaiterObject(eventObject);}
void WaitForEvent(EventObject *eventObject) {WaitForSingleWaiterObject(eventObject);}
bool WaitForEventWithTimeout(EventObject *eventObject, _int64 timeoutInMillis)
{
    DWORD retVal = WaitForSingleObjectEx(*eventObject, (unsigned)timeoutInMillis, FALSE);

    if (retVal == WAIT_OBJECT_0) {
        return true;
    } else if (retVal == WAIT_TIMEOUT) {
        return false;
    }
    WriteErrorMessage("WaitForSingleObject returned unexpected result %d (error %d)\n", retVal, GetLastError());
    soft_exit(1);
    return false;   // NOTREACHED: Just to avoid the compiler complaining.
}


void BindThreadToProcessor(unsigned processorNumber) // This hard binds a thread to a processor.  You can no-op it at some perf hit.
{
    if (!SetThreadAffinityMask(GetCurrentThread(),((unsigned _int64)1) << processorNumber)) {
        WriteErrorMessage("Binding thread to processor %d failed, %d\n",processorNumber,GetLastError());
    }
}

int InterlockedIncrementAndReturnNewValue(volatile int *valueToIncrement)
{
    return InterlockedIncrement((volatile long *)valueToIncrement);
}

int InterlockedDecrementAndReturnNewValue(volatile int *valueToDecrement)
{
    return InterlockedDecrement((volatile long *)valueToDecrement);
}

_uint32 InterlockedCompareExchange32AndReturnOldValue(volatile _uint32 *valueToUpdate, _uint32 replacementValue, _uint32 desiredPreviousValue)
{
    return (_uint32) InterlockedCompareExchange(valueToUpdate, replacementValue, desiredPreviousValue);
}

_uint64 InterlockedCompareExchange64AndReturnOldValue(volatile _uint64 *valueToUpdate, _uint64 replacementValue, _uint64 desiredPreviousValue)
{
    return (_uint64) InterlockedCompareExchange(valueToUpdate, replacementValue, desiredPreviousValue);
}

void* InterlockedCompareExchangePointerAndReturnOldValue(void * volatile *valueToUpdate, void* replacementValue, void* desiredPreviousValue)
{
    return InterlockedCompareExchangePointer(valueToUpdate, replacementValue, desiredPreviousValue);
}

struct WrapperThreadContext {
    ThreadMainFunction      mainFunction;
    void                    *mainFunctionParameter;
};

    DWORD WINAPI
WrapperThreadMain(PVOID Context)
{
    WrapperThreadContext *context = (WrapperThreadContext *)Context;

    (*context->mainFunction)(context->mainFunctionParameter);
    delete context;
    context = NULL;
    return 0;
}


bool StartNewThread(ThreadMainFunction threadMainFunction, void *threadMainFunctionParameter)
{
    WrapperThreadContext *context = new WrapperThreadContext;
    if (NULL == context) {
        return false;
    }
    context->mainFunction = threadMainFunction;
    context->mainFunctionParameter = threadMainFunctionParameter;

    HANDLE hThread;
    DWORD threadId;
    hThread = CreateThread(NULL,0,WrapperThreadMain,context,0,&threadId);
    if (NULL == hThread) {
        WriteErrorMessage("Create thread failed, %d\n",GetLastError());
        delete context;
        context = NULL;
        return false;
    }

    CloseHandle(hThread);
    hThread = NULL;
    return true;
}

void SleepForMillis(unsigned millis)
{
  Sleep(millis);
}

unsigned GetNumberOfProcessors()
{
    SYSTEM_INFO systemInfo[1];
    GetSystemInfo(systemInfo);

    return systemInfo->dwNumberOfProcessors;
}

_int64 QueryFileSize(const char *fileName) {
    HANDLE hFile = CreateFile(fileName,GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,NULL);
    if (INVALID_HANDLE_VALUE == hFile) {
        WriteErrorMessage("Unable to open file '%s' for QueryFileSize, %d\n", fileName, GetLastError());
        soft_exit(1);
    }
    LARGE_INTEGER fileSize;
    if (!GetFileSizeEx(hFile,&fileSize)) {
        WriteErrorMessage("GetFileSize failed, %d\n",GetLastError());
        soft_exit(1);
    }
    CloseHandle(hFile);

    return fileSize.QuadPart;
}

    bool
DeleteSingleFile(
    const char* filename)
{
    return DeleteFile(filename) ? true : false;
}

    bool
MoveSingleFile(
    const char* oldFileName,
    const char* newFileName)
{
    return MoveFile(oldFileName, newFileName) ? true : false;
}

class LargeFileHandle
{
public:
    HANDLE  handle;
};

    LargeFileHandle*
OpenLargeFile(
    const char* filename,
    const char* mode)
{
    _ASSERT(strlen(mode) == 1 && (*mode == 'r' || *mode == 'w' || *mode == 'a'));
    LargeFileHandle* result = new LargeFileHandle();
    result->handle = CreateFile(filename,
        *mode == 'r' ? GENERIC_READ :
            *mode == 'a' ? FILE_APPEND_DATA
            : GENERIC_WRITE,
        0 /* exclusive */,
        NULL,
        *mode == 'w' ? CREATE_ALWAYS : OPEN_EXISTING,
        FILE_FLAG_SEQUENTIAL_SCAN,
        NULL);
    if (result->handle == NULL) {
        WriteErrorMessage("open large file %s failed with 0x%x\n", filename, GetLastError());
        delete (void*) result;
        return NULL;
    }
    return result;
}

    size_t
WriteLargeFile(
    LargeFileHandle* file,
    void* buffer,
    size_t bytes)
{
    size_t count = bytes;
    while (count > 0) {
        DWORD step = 0;
        if ((! WriteFile(file->handle, buffer, (DWORD) min(count, (size_t) 0x2000000), &step, NULL)) || step == 0) {
            WriteErrorMessage("WriteLargeFile failed at %lu of %lu bytes with 0x%x\n", bytes - count, bytes, GetLastError());
            return bytes - count;
        }
        count -= step;
        buffer = ((char*) buffer) + step;
    }
    return bytes;
}


    size_t
ReadLargeFile(
    LargeFileHandle* file,
    void* buffer,
    size_t bytes)
{
    size_t count = bytes;
    while (count > 0) {
        DWORD step = 0;
        if ((! ReadFile(file->handle, buffer, (DWORD) min(count, (size_t) 0x1000000), &step, NULL)) || step == 0) {
            WriteErrorMessage("ReadLargeFile failed at %lu of %lu bytes with 0x%x\n", bytes - count, bytes, GetLastError());
            return bytes - count;
        }
        count -= step;
        buffer = ((char*) buffer) + step;
    }
    return bytes;
}

    void
CloseLargeFile(
    LargeFileHandle* file)
{
    if (CloseHandle(file->handle)) {
        delete (void*) file;
    }
}

class MemoryMappedFile
{
public:
    HANDLE  fileHandle;
    HANDLE  fileMapping;
    void*   mappedAddress;
};


    MemoryMappedFile*
OpenMemoryMappedFile(
    const char* filename,
    size_t offset,
    size_t length,
    void** o_contents,
    bool write,
    bool sequential)
{
    MemoryMappedFile* result = new MemoryMappedFile();
    result->fileHandle = CreateFile(filename, (write ? GENERIC_WRITE : 0) | GENERIC_READ, 0, NULL, OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL | (sequential ? FILE_FLAG_SEQUENTIAL_SCAN : FILE_FLAG_RANDOM_ACCESS), NULL);
    if (result->fileHandle == NULL) {
        WriteErrorMessage("unable to open mapped file %s error 0x%x\n", filename, GetLastError());
        delete result;
        return NULL;
    }
    result->fileMapping = CreateFileMapping(result->fileHandle, NULL, write ? PAGE_READWRITE : PAGE_READONLY, 0, 0, NULL);
    if (result->fileMapping == NULL) {
        WriteErrorMessage("unable to create file mapping %s error 0x%x\n", filename, GetLastError());
        delete result;
        return NULL;
    }
    *o_contents = result->mappedAddress = MapViewOfFile(result->fileMapping,
        write ? FILE_MAP_WRITE : FILE_MAP_READ,
        (DWORD) (offset >> (8 * sizeof(DWORD))),
        (DWORD) offset,
        length);
    if (*o_contents == NULL) {
        WriteErrorMessage("unable to map file %s error 0x%x\n", filename, GetLastError());
        delete result;
        return NULL;
    }
    return result;
}

    void
CloseMemoryMappedFile(
    MemoryMappedFile* mappedFile)
{
    bool ok = UnmapViewOfFile(mappedFile->mappedAddress) &&
        CloseHandle(mappedFile->fileMapping) &&
        CloseHandle(mappedFile->fileHandle);
    if (ok) {
        delete (void*) mappedFile;
    } else {
        WriteErrorMessage("unable to close memory mapped file, error 0x%x\n", GetLastError());
    }
}

void AdviseMemoryMappedFilePrefetch(const MemoryMappedFile *mappedFile)
{
  // No-op on WIndows.
}


class WindowsAsyncFile : public AsyncFile
{
public:
    static WindowsAsyncFile* open(const char* filename, bool write);

    WindowsAsyncFile(HANDLE i_hFile);

    virtual bool close();

    class Writer : public AsyncFile::Writer
    {
    public:
        Writer(WindowsAsyncFile* i_file);

        virtual bool close();

        virtual bool beginWrite(void* buffer, size_t length, size_t offset, size_t *bytesWritten);

        virtual bool waitForCompletion();
    
    private:
        WindowsAsyncFile*   file;
        bool                writing;
        OVERLAPPED          lap;
    };

    virtual AsyncFile::Writer* getWriter();
    
    class Reader : public AsyncFile::Reader
    {
    public:
        Reader(WindowsAsyncFile* i_file);

        virtual bool close();

        virtual bool beginRead(void* buffer, size_t length, size_t offset, size_t *bytesRead);

        virtual bool waitForCompletion();
    
    private:
        WindowsAsyncFile*   file;
        bool                reading;
        OVERLAPPED          lap;
    };

    virtual AsyncFile::Reader* getReader();

private:
    HANDLE      hFile;
};

    WindowsAsyncFile*
WindowsAsyncFile::open(
    const char* filename,
    bool write)
{
    HANDLE hFile = CreateFile(filename,
        GENERIC_READ | (write ? GENERIC_WRITE : 0),
        write ? 0 : FILE_SHARE_READ,
        NULL,
        write ? CREATE_ALWAYS : OPEN_EXISTING,
        FILE_FLAG_OVERLAPPED,
        NULL);
    if (INVALID_HANDLE_VALUE == hFile) {
        WriteErrorMessage("Unable to create SAM file '%s', %d\n",filename,GetLastError());
        return NULL;
    }
    return new WindowsAsyncFile(hFile);
}

WindowsAsyncFile::WindowsAsyncFile(
    HANDLE i_hFile)
    : hFile(i_hFile)
{
}

    bool
WindowsAsyncFile::close()
{
    return CloseHandle(hFile) ? true : false;
}

    AsyncFile::Writer*
WindowsAsyncFile::getWriter()
{
    return new Writer(this);
}

WindowsAsyncFile::Writer::Writer(WindowsAsyncFile* i_file)
    : file(i_file), writing(false)
{
    lap.hEvent = CreateEvent(NULL,FALSE,FALSE,NULL);
}

    bool
WindowsAsyncFile::Writer::close()
{
    waitForCompletion();
    return CloseHandle(lap.hEvent) ? true : false;
}

    bool
WindowsAsyncFile::Writer::beginWrite(
    void* buffer,
    size_t length,
    size_t offset,
    size_t *bytesWritten)
{
    if (! waitForCompletion()) {
        return false;
    }
    lap.OffsetHigh = (DWORD) (offset >> (8 * sizeof(DWORD)));
    lap.Offset = (DWORD) offset;
    if (!WriteFile(file->hFile,buffer, (DWORD) length, (LPDWORD) bytesWritten, &lap)) {
        if (ERROR_IO_PENDING != GetLastError()) {
            WriteErrorMessage("WindowsAsyncFile: WriteFile failed, %d\n",GetLastError());
            return false;
        }
    }
    writing = true;
    return true;
}

    bool
WindowsAsyncFile::Writer::waitForCompletion()
{
    if (writing) {
        DWORD nBytesTransferred;
        if (!GetOverlappedResult(file->hFile,&lap,&nBytesTransferred,TRUE)) {
            return false;
        }
        writing = false;
    }
    return true;
}

    AsyncFile::Reader*
WindowsAsyncFile::getReader()
{
    return new Reader(this);
}

WindowsAsyncFile::Reader::Reader(
    WindowsAsyncFile* i_file)
    : file(i_file), reading(false)
{
    lap.hEvent = CreateEvent(NULL,FALSE,FALSE,NULL);
}

    bool
WindowsAsyncFile::Reader::close()
{
    return CloseHandle(lap.hEvent) ? true : false;
}

    bool
WindowsAsyncFile::Reader::beginRead(
    void* buffer,
    size_t length,
    size_t offset,
    size_t* bytesRead)
{
    if (! waitForCompletion()) {
        return false;
    }
    lap.OffsetHigh = (DWORD) (offset >> (8 * sizeof(DWORD)));
    lap.Offset = (DWORD) offset;
    if (!ReadFile(file->hFile, buffer,(DWORD) length, (LPDWORD) bytesRead, &lap)) {
        if (ERROR_IO_PENDING != GetLastError()) {
            WriteErrorMessage("WindowsSAMWriter: WriteFile failed, %d\n",GetLastError());
            return false;
        }
    }
    reading = true;
    return true;
}

    bool
WindowsAsyncFile::Reader::waitForCompletion()
{
    if (reading) {
        DWORD nBytesTransferred;
        if (!GetOverlappedResult(file->hFile,&lap,&nBytesTransferred,TRUE)) {
            return false;
        }
        reading = false;
    }
    return true;
}

_int64 InterlockedAdd64AndReturnNewValue(volatile _int64 *valueToWhichToAdd, _int64 amountToAdd)
{
    return InterlockedAdd64((volatile LONGLONG *)valueToWhichToAdd,(LONGLONG)amountToAdd);
}

int _fseek64bit(FILE *stream, _int64 offset, int origin)
{
    return _fseeki64(stream,offset,origin);
}

int getpagesize()
{
    SYSTEM_INFO systemInfo;
    GetSystemInfo(&systemInfo);
    return systemInfo.dwAllocationGranularity;
}

FileMapper::FileMapper()
{
    hFile = INVALID_HANDLE_VALUE;
    hMapping = NULL;
    initialized = false;
    pagesize = getpagesize();
    mapCount = 0;
    
#if 0
    hFilePrefetch = INVALID_HANDLE_VALUE;
    lap->hEvent = NULL;
    prefetchBuffer = BigAlloc(prefetchBufferSize);
    isPrefetchOutstanding = false;
    lastPrefetch = 0;
#endif

    millisSpentInReadFile = 0;
    countOfImmediateCompletions = 0;
    countOfDelayedCompletions = 0;
    countOfFailures = 0;
}

bool
FileMapper::init(const char *i_fileName)
{
    if (initialized) {
        if (strcmp(fileName, i_fileName)) {
            WriteErrorMessage("FileMapper already initialized with %s, cannot init with %s\n", fileName, i_fileName);
            return false;
        }
        return true;
    }
    fileName = i_fileName;
    hFile = CreateFile(fileName, GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,FILE_FLAG_SEQUENTIAL_SCAN,NULL);
    if (INVALID_HANDLE_VALUE == hFile) {
        WriteErrorMessage("Failed to open '%s', error %d\n",fileName, GetLastError());
        return false;
    }

#if 0
    hFilePrefetch = CreateFile(fileName, GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,FILE_FLAG_OVERLAPPED,NULL);
    if (INVALID_HANDLE_VALUE == hFilePrefetch) {
        WriteErrorMessage("Failed to open '%s' for prefetch, error %d\n",fileName, GetLastError());
        CloseHandle(hFile);
        return false;
    }
#endif

    BY_HANDLE_FILE_INFORMATION fileInfo;
    if (!GetFileInformationByHandle(hFile,&fileInfo)) {
        WriteErrorMessage("Unable to get file information for '%s', error %d\n", fileName, GetLastError());
        CloseHandle(hFile);
#if 0
        CloseHandle(hFilePrefetch);
#endif
        return false;
    }
    LARGE_INTEGER liFileSize;
    liFileSize.HighPart = fileInfo.nFileSizeHigh;
    liFileSize.LowPart = fileInfo.nFileSizeLow;
    fileSize = liFileSize.QuadPart;

    hMapping = CreateFileMapping(hFile,NULL,PAGE_READONLY,0,0,NULL);
    if (NULL == hMapping) {
        WriteErrorMessage("Unable to create mapping to file '%s', %d\n", fileName, GetLastError());
        CloseHandle(hFile);
#if 0
        CloseHandle(hFilePrefetch);
#endif
        return false;
    }
#if 0
    lap->hEvent = CreateEvent(NULL,FALSE,FALSE,NULL);
#endif
    initialized = true;

 
    return true;
}

char *
FileMapper::createMapping(size_t offset, size_t amountToMap, void** o_mappedBase)
{
    size_t beginRounding = offset % pagesize;
    LARGE_INTEGER liStartingOffset;
    liStartingOffset.QuadPart = offset - beginRounding;

    size_t endRounding = 0;
    if ((amountToMap + beginRounding) % pagesize != 0) {
        endRounding = pagesize - (amountToMap + beginRounding) % pagesize;
    }
    size_t mapRequestSize = beginRounding + amountToMap + endRounding;
    _ASSERT(mapRequestSize % pagesize == 0);
    if (mapRequestSize + liStartingOffset.QuadPart >= fileSize) {
        mapRequestSize = 0; // Says to just map the whole thing.
    }

    char* mappedBase = (char *)MapViewOfFile(hMapping,FILE_MAP_READ,liStartingOffset.HighPart,liStartingOffset.LowPart, mapRequestSize);
    if (NULL == mappedBase) {
        WriteErrorMessage("Unable to map file, %d\n", GetLastError());
        return NULL;
    } 
    char* mappedRegion = mappedBase + beginRounding;

#if 0
    prefetch(0);
#endif

    InterlockedIncrementAndReturnNewValue(&mapCount);
    *o_mappedBase = mappedBase;
   return mappedRegion;
}

void
FileMapper::unmap(void* mappedBase)
{
    _ASSERT(mapCount > 0);
    if (mapCount > 0) {
        int n = InterlockedDecrementAndReturnNewValue(&mapCount);
        _ASSERT(n >= 0);
        if (!UnmapViewOfFile(mappedBase)) {
            WriteErrorMessage("Unmap of file failed, %d\n", GetLastError());
        }
    }
}

FileMapper::~FileMapper()
{
    _ASSERT(mapCount == 0);
#if 0
    if (isPrefetchOutstanding) {
        DWORD numberOfBytesTransferred;
        GetOverlappedResult(hFile,lap,&numberOfBytesTransferred,TRUE);
    }
    BigDealloc(prefetchBuffer);
    prefetchBuffer = NULL;
    CloseHandle(hFilePrefetch);
    CloseHandle(lap->hEvent);
#endif
    CloseHandle(hMapping);
    CloseHandle(hFile);
    WriteErrorMessage("FileMapper: %lld immediate completions, %lld delayed completions, %lld failures, %lld ms in readfile (%lld ms/call)\n",countOfImmediateCompletions, countOfDelayedCompletions, countOfFailures, millisSpentInReadFile, 
        millisSpentInReadFile/max(1, countOfImmediateCompletions + countOfDelayedCompletions + countOfFailures));
}
    
#if 0
void
FileMapper::prefetch(size_t currentRead)
{
    if (currentRead + prefetchBufferSize / 2 <= lastPrefetch || lastPrefetch + prefetchBufferSize >= amountMapped) {
        //
        // Nothing to do; we're either not ready for more prefetching or we're at the end of our region.
        //
        return;
    }

    if (isPrefetchOutstanding) {
        //
        // See if the last prefetch is done.
        //
        DWORD numberOfBytesTransferred;
        if (GetOverlappedResult(hFile,lap,&numberOfBytesTransferred,FALSE)) {
            isPrefetchOutstanding = false;
        } else {
#if     DBG
            if (GetLastError() != ERROR_IO_PENDING) {
                WriteErrorMessage("mapped file prefetcher: GetOverlappedResult failed, %d\n", GetLastError());
            }
#endif  // DBG
            return;  // There's still IO on outstanding, we can't start more.
        }
    }

    DWORD amountToRead = (DWORD)__min(prefetchBufferSize, amountMapped - lastPrefetch);
    _ASSERT(amountToRead > 0);  // Else we should have failed the initial check and returned
    LARGE_INTEGER liReadOffset;
    lastPrefetch += prefetchBufferSize;
    liReadOffset.QuadPart = lastPrefetch;
    lap->OffsetHigh = liReadOffset.HighPart;
    lap->Offset = liReadOffset.LowPart;
    DWORD nBytesRead;

    _int64 start = timeInMillis();
    if (!ReadFile(hFilePrefetch,prefetchBuffer,amountToRead,&nBytesRead,lap)) {
        if (GetLastError() == ERROR_IO_PENDING) {
            InterlockedAdd64AndReturnNewValue(&countOfDelayedCompletions,1);
            isPrefetchOutstanding = true;
        } else {
           InterlockedAdd64AndReturnNewValue(&countOfFailures,1);
#if     DBG
            if (GetLastError() != ERROR_IO_PENDING) {
                WriteErrorMessage("mapped file prefetcher: ReadFile failed, %d\n", GetLastError());
            }
#endif  // DBG
            isPrefetchOutstanding = false; // Just ignore it
        }
    } else {
        InterlockedAdd64AndReturnNewValue(&countOfImmediateCompletions,1);
        isPrefetchOutstanding = false;
    }
    InterlockedAdd64AndReturnNewValue(&millisSpentInReadFile,timeInMillis() - start);
}
#endif

void PreventMachineHibernationWhileThisThreadIsAlive()
{
	SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED);
}

void SetToLowSchedulingPriority()
{
    if (!SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS)) {
        WriteErrorMessage("Unable to set process to background priority class, %d.  Ignoring and proceeding at normal priority\n", GetLastError());
    }
}

struct NamedPipe {
	HANDLE hPipe;
};

NamedPipe *OpenNamedPipe(const char *pipeName, bool serverSide)
{
	NamedPipe *pipe = new NamedPipe;
	const char *prefix = "\\\\.\\pipe\\";

	char *fullyQualifiedPipeName = new char[strlen(prefix) + strlen(pipeName) + 1];	// +1 for null

	sprintf(fullyQualifiedPipeName, "%s%s", prefix, pipeName);

	if (serverSide) {
		pipe->hPipe = CreateNamedPipe(fullyQualifiedPipeName, PIPE_ACCESS_DUPLEX | FILE_FLAG_FIRST_PIPE_INSTANCE, PIPE_TYPE_MESSAGE | PIPE_READMODE_MESSAGE | PIPE_WAIT | PIPE_REJECT_REMOTE_CLIENTS, 1, 100000, 100000, 0, NULL);
		if (INVALID_HANDLE_VALUE == pipe->hPipe) {
			WriteErrorMessage("OpenNamedPipe('%s'): unable to open pipe, %d\n", fullyQualifiedPipeName, GetLastError());
			delete pipe;
			return NULL;
		}

		if (!ConnectNamedPipe(pipe->hPipe, NULL) && ERROR_PIPE_CONNECTED != GetLastError()) {
			WriteErrorMessage("Unable to connect named pipe, %d\n", GetLastError());
			delete pipe;
			return NULL;
		}
	} else {
		pipe->hPipe = CreateFile(fullyQualifiedPipeName, GENERIC_READ | GENERIC_WRITE, 0, NULL, OPEN_EXISTING, 0, NULL);
		while (INVALID_HANDLE_VALUE == pipe->hPipe && GetLastError() == ERROR_PIPE_BUSY) {
			WriteStatusMessage("Server is busy.  Waiting.\n");
			if (!WaitNamedPipe(fullyQualifiedPipeName, NMPWAIT_WAIT_FOREVER)) {
				fprintf(stderr, "Waiting for server connection failed, %d\n", GetLastError());
				delete pipe;
				return NULL;
			}
			pipe->hPipe = CreateFile(fullyQualifiedPipeName, GENERIC_READ | GENERIC_WRITE, 0, NULL, OPEN_EXISTING, 0, NULL);
		}

		if (INVALID_HANDLE_VALUE == pipe->hPipe) {
			WriteErrorMessage("Unable to open server connection '%s', %d\n", fullyQualifiedPipeName, GetLastError());
			delete pipe;
			return NULL;
		}
	}

	return pipe;
}

bool ReadFromNamedPipe(NamedPipe *pipe, char *outputBuffer, size_t outputBufferSize)
{
	DWORD bytesRead;
	while (!ReadFile(pipe->hPipe, outputBuffer, (DWORD)outputBufferSize, &bytesRead, NULL)) {
		if (GetLastError() != ERROR_BROKEN_PIPE && GetLastError() != ERROR_NO_DATA) {
			fprintf(stderr, "Read named pipe failed, %d\n", GetLastError());	// Don't use WriteErrorMessage, it will try to send on the pipe
			return false;
		}

		if (GetLastError() == ERROR_BROKEN_PIPE) {
			if (!DisconnectNamedPipe(pipe->hPipe)) {
				fprintf(stderr, "Disconnect named pipe failed, %d; ignoring\n", GetLastError());
			}
			if (!ConnectNamedPipe(pipe->hPipe, NULL) && ERROR_PIPE_CONNECTED != GetLastError()) {
				fprintf(stderr, "ReadFromNamedPipe: reconnecting to pipe failed, %d\n", GetLastError());
				return false;
			}
		}
	}
	return true;
}

bool WriteToNamedPipe(NamedPipe *pipe, const char *stringToWrite)
{
	DWORD bytesWritten;
	if (!WriteFile(pipe->hPipe, stringToWrite, (DWORD)strlen(stringToWrite) + 1, &bytesWritten, NULL)) {	// +1 sends terminating NULL
		fprintf(stderr, "WriteToNamedPipe: write failed, %d\n", GetLastError());
		return false;
	}

	FlushFileBuffers(pipe->hPipe);

	if (bytesWritten != strlen(stringToWrite) + 1) {
		fprintf(stderr, "WriteToNamedPipe:  expected to write %lld bytes, actually wrote %d\n", strlen(stringToWrite) + 1, bytesWritten);
	}

	return bytesWritten == strlen(stringToWrite) + 1;
}

void CloseNamedPipe(NamedPipe *pipe)
{
	CloseHandle(pipe->hPipe);
	delete pipe;
}

const char *DEFAULT_NAMED_PIPE_NAME = "SNAP";
#else   // _MSC_VER

#if defined(__MACH__)
#include <mach/clock.h>
#include <mach/mach.h>
#endif

_int64 timeInMillis()
/**
 * Get the current time in milliseconds since some arbitrary starting point
 * (e.g. system boot or epoch)
 */
{
    timeval t;
    gettimeofday(&t, NULL);
    return ((_int64) t.tv_sec) * 1000 + ((_int64) t.tv_usec) / 1000;
}

_int64 timeInNanos()
{
    timespec ts;
#if defined(__linux__)
    clock_gettime(CLOCK_REALTIME, &ts); // Works on Linux
#elif defined(__MACH__)
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts.tv_sec = mts.tv_sec;
    ts.tv_nsec = mts.tv_nsec;
#else
    #error "Don't know how to get time in nanos on your platform"
#endif
    return ((_int64) ts.tv_sec) * 1000000000 + (_int64) ts.tv_nsec;
}

void AcquireUnderlyingExclusiveLock(UnderlyingExclusiveLock *lock)
{
    pthread_mutex_lock(lock);
}

void ReleaseUnderlyingExclusiveLock(UnderlyingExclusiveLock *lock)
{
    pthread_mutex_unlock(lock);
}

bool InitializeUnderlyingExclusiveLock(UnderlyingExclusiveLock *lock)
{
    return pthread_mutex_init(lock, NULL) == 0;
}

bool DestroyUnderlyingExclusiveLock(UnderlyingExclusiveLock *lock)
{
    return pthread_mutex_destroy(lock) == 0;
}

class SingleWaiterObjectImpl {
protected:
    pthread_mutex_t lock;
    pthread_cond_t cond;
    bool set;

public:
    bool init() {
        if (pthread_mutex_init(&lock, NULL) != 0) {
            return false;
        }
        if (pthread_cond_init(&cond, NULL) != 0) {
            pthread_mutex_destroy(&lock);
            return false;
        }
        set = false;
        return true;
    }

    void signal() {
        pthread_mutex_lock(&lock);
        set = true;
        pthread_cond_signal(&cond);
        pthread_mutex_unlock(&lock);
    }

    void wait() {
        pthread_mutex_lock(&lock);
        while (!set) {
            pthread_cond_wait(&cond, &lock);
        }
        pthread_mutex_unlock(&lock);
    }

    bool waitWithTimeout(_int64 timeoutInMillis) {
        struct timespec wakeTime;
#ifdef __LINUX__
        clock_gettime(CLOCK_REALTIME, &wakeTime);
        wakeTime.tv_nsec += timeoutInMillis * 1000000;
#elif defined(__MACH__)
        clock_serv_t cclock;
        mach_timespec_t mts;
        host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
        clock_get_time(cclock, &mts);
        mach_port_deallocate(mach_task_self(), cclock);
        wakeTime.tv_nsec = mts.tv_nsec + timeoutInMillis * 1000000;
        wakeTime.tv_sec = mts.tv_sec;
#endif
        wakeTime.tv_sec += wakeTime.tv_nsec / 1000000000;
        wakeTime.tv_nsec = wakeTime.tv_nsec % 1000000000;

        bool timedOut = false;

        pthread_mutex_lock(&lock);
        while (!set) {
            int retVal = pthread_cond_timedwait(&cond, &lock, &wakeTime);
            if (retVal == ETIMEDOUT) {
                timedOut = true;
                break;
            }
        }
        pthread_mutex_unlock(&lock);
        return !timedOut;
    }

    bool destroy() {
        pthread_cond_destroy(&cond);
        pthread_mutex_destroy(&lock);
    }
};

bool CreateSingleWaiterObject(SingleWaiterObject *waiter)
{
    SingleWaiterObjectImpl *obj = new SingleWaiterObjectImpl;
    if (obj == NULL) {
        return false;
    }
    if (!obj->init()) {
        delete obj;
        return false;
    }
    *waiter = obj;
    return true;
}

void DestroySingleWaiterObject(SingleWaiterObject *waiter)
{
    (*waiter)->destroy();
    delete *waiter;
}

void SignalSingleWaiterObject(SingleWaiterObject *waiter)
{
    (*waiter)->signal();
}

bool WaitForSingleWaiterObject(SingleWaiterObject *waiter)
{
    (*waiter)->wait();
    return true;
}

void ResetSingleWaiterObject(SingleWaiterObject *waiter)
{
    (*waiter)->init();
}

class EventObjectImpl : public SingleWaiterObjectImpl
{
public:
    void signalAll()
    {
        pthread_mutex_lock(&lock);
        set = true;
        pthread_cond_broadcast(&cond);
        pthread_mutex_unlock(&lock);
    }
    void blockAll()
    {
        pthread_mutex_lock(&lock);
	    set = false;
        pthread_mutex_unlock(&lock);
    }
};

void CreateEventObject(EventObject *newEvent)
{
    EventObjectImpl* obj = new EventObjectImpl();
    if (obj == NULL) {
        return;
    }
    if (!obj->init()) {
        delete obj;
        return;
    }
    *newEvent = obj;
}

void DestroyEventObject(EventObject *eventObject)
{
    (*eventObject)->destroy();
    delete *eventObject;
}

void AllowEventWaitersToProceed(EventObject *eventObject)
{
    (*eventObject)->signalAll();
}

void PreventEventWaitersFromProceeding(EventObject *eventObject)
{
  (*eventObject)->blockAll();
}

void WaitForEvent(EventObject *eventObject)
{
  (*eventObject)->wait();
}

bool WaitForEventWithTimeout(EventObject *eventObject, _int64 timeoutInMillis)
{
    return (*eventObject)->waitWithTimeout(timeoutInMillis);
}

int InterlockedIncrementAndReturnNewValue(volatile int *valueToDecrement)
{
    return (int) __sync_fetch_and_add((volatile int*) valueToDecrement, 1) + 1;
}

int InterlockedDecrementAndReturnNewValue(volatile int *valueToDecrement)
{
    return __sync_fetch_and_sub(valueToDecrement, 1) - 1;
}

_int64 InterlockedAdd64AndReturnNewValue(volatile _int64 *valueToWhichToAdd, _int64 amountToAdd)
{
    return __sync_fetch_and_add(valueToWhichToAdd, amountToAdd) + amountToAdd;
}

_uint32 InterlockedCompareExchange32AndReturnOldValue(volatile _uint32 *valueToUpdate, _uint32 replacementValue, _uint32 desiredPreviousValue)
{
  return __sync_val_compare_and_swap(valueToUpdate, desiredPreviousValue, replacementValue);
}

_uint64 InterlockedCompareExchange64AndReturnOldValue(volatile _uint64 *valueToUpdate, _uint64 replacementValue, _uint64 desiredPreviousValue)
{
  return (_uint64) __sync_val_compare_and_swap((volatile _int64 *) valueToUpdate, desiredPreviousValue, replacementValue);
}

void* InterlockedCompareExchangePointerAndReturnOldValue(void * volatile *valueToUpdate, void* replacementValue, void* desiredPreviousValue)
{
  return __sync_val_compare_and_swap(valueToUpdate, desiredPreviousValue, replacementValue);
}

namespace {

// POSIX thread functions need to return void*, so we wrap the ThreadMainFunction in our API
struct ThreadInfo {
    ThreadMainFunction function;
    void *parameter;

    ThreadInfo(ThreadMainFunction f, void *p): function(f), parameter(p) {}
};

void* runThread(void* infoVoidPtr) {
    ThreadInfo *info = (ThreadInfo*) infoVoidPtr;
    info->function(info->parameter);
    delete info;
    return NULL;
}

}

bool StartNewThread(ThreadMainFunction threadMainFunction, void *threadMainFunctionParameter)
{
    ThreadInfo *info = new ThreadInfo(threadMainFunction, threadMainFunctionParameter);
    pthread_t thread;
    return pthread_create(&thread, NULL, runThread, info) == 0;
}

void BindThreadToProcessor(unsigned processorNumber)
{
#ifdef __linux__
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(processorNumber, &cpuset);
    if (sched_setaffinity(0, sizeof(cpu_set_t), &cpuset) != 0) {
        perror("sched_setaffinity");
    }
#endif
}

unsigned GetNumberOfProcessors()
{
    return (unsigned) sysconf(_SC_NPROCESSORS_ONLN);
}

void SleepForMillis(unsigned millis)
{
  usleep(millis*1000);
}

_int64 QueryFileSize(const char *fileName)
{
    int fd = open(fileName, O_RDONLY);
    _ASSERT(fd != -1);
    struct stat sb;
    int r = fstat(fd, &sb);
    _ASSERT(r != -1);
    _int64 fileSize = sb.st_size;
    close(fd);
    return fileSize;
}

    bool
DeleteSingleFile(
    const char* filename)
{
    return unlink(filename) == 0;
}

    bool
MoveSingleFile(
    const char* from,
    const char* to)
{
    return rename(from, to) == 0;
}

class LargeFileHandle
{
public:
    FILE* file;
};

    LargeFileHandle*
OpenLargeFile(
    const char* filename,
    const char* mode)
{
    _ASSERT(strlen(mode) == 1 && (*mode == 'r' || *mode == 'w' || *mode == 'a'));
    char fmode[3];
    fmode[0] = mode[0]; fmode[1] = 'b'; fmode[2] = '\0';
    FILE* file = fopen(filename, fmode);
    if (file == NULL) {
        return NULL;
    }
    LargeFileHandle* result = new LargeFileHandle();
    result->file = file;
    return result;
}

    size_t
WriteLargeFile(
    LargeFileHandle* file,
    void* buffer,
    size_t bytes)
{
    return fwrite(buffer, 1, bytes, file->file);
}


    size_t
ReadLargeFile(
    LargeFileHandle* file,
    void* buffer,
    size_t bytes)
{
    return fread(buffer, 1, bytes, file->file);    
}

    void
CloseLargeFile(
    LargeFileHandle* file)
{
    if (0 == fclose(file->file)) {
        delete file;
    }
}

class MemoryMappedFile
{
public:
    int     fd;
    void*   map;
    size_t  length;
};

    MemoryMappedFile*
OpenMemoryMappedFile(
    const char* filename,
    size_t offset,
    size_t length,
    void** o_contents,
    bool write,
    bool sequential)
{
  int fd = open(filename, write ? O_CREAT | O_RDWR : O_RDONLY, S_IRUSR | S_IWUSR);
    if (fd < 0) {
        warn("OpenMemoryMappedFile %s failed", filename);
        return NULL;
    }
    // todo: large page support
    size_t page = getpagesize();
    size_t extra = offset % page;
    void* map = mmap(NULL, length + extra, (write ? PROT_WRITE : 0) | PROT_READ, MAP_PRIVATE, fd, offset - extra);
    if (map == NULL || map == MAP_FAILED) {
        warn("OpenMemoryMappedFile %s mmap failed", filename);
        close(fd);
        return NULL;
    }
    int e = madvise(map, length + extra, sequential ? MADV_SEQUENTIAL : MADV_RANDOM);
    if (e < 0) {
        warn("OpenMemoryMappedFile %s madvise failed", filename);
    }
    MemoryMappedFile* result = new MemoryMappedFile();
    result->fd = fd;
    result->map = map;
    result->length = length + extra;
    *o_contents = (char*)map + extra;
    return result;
}

    void
CloseMemoryMappedFile(
    MemoryMappedFile* mappedFile)
{
    int e = munmap(mappedFile->map, mappedFile->length);
    int e2 = close(mappedFile->fd);
    if (e != 0 || e2 != 0) {
        WriteErrorMessage("CloseMemoryMapped file failed\n");
    }
}

void AdviseMemoryMappedFilePrefetch(const MemoryMappedFile *mappedFile)
{
  if (madvise(mappedFile->map, mappedFile->length, MADV_SEQUENTIAL)) {
    WriteErrorMessage("madvise MADV_SEQUENTIAL failed (since it's only an optimization, this is OK).  Errno %d\n", errno);
  }

  if (madvise(mappedFile->map, mappedFile->length, MADV_WILLNEED)) {
    WriteErrorMessage("madvise MADV_WILLNEED failed (since it's only an optimization, this is OK).  Errno %d\n", errno);
  }
}

#ifdef __linux__

class PosixAsyncFile : public AsyncFile
{
public:
    static PosixAsyncFile* open(const char* filename, bool write);

    PosixAsyncFile(int i_fd);

    virtual bool close();

    class Writer : public AsyncFile::Writer
    {
    public:
        Writer(PosixAsyncFile* i_file);

        virtual bool close();

        virtual bool beginWrite(void* buffer, size_t length, size_t offset, size_t *bytesWritten);

        virtual bool waitForCompletion();
    
    private:
        PosixAsyncFile*     file;
        bool                writing;
        SingleWaiterObject  ready;
        struct aiocb        aiocb;
        size_t*             result;
    };

    virtual AsyncFile::Writer* getWriter();
    
    class Reader : public AsyncFile::Reader
    {
    public:
        Reader(PosixAsyncFile* i_file);

        virtual bool close();

        virtual bool beginRead(void* buffer, size_t length, size_t offset, size_t *bytesRead);

        virtual bool waitForCompletion();
    
    private:
        PosixAsyncFile*     file;
        bool                reading;
        SingleWaiterObject  ready;
        struct aiocb        aiocb;
        size_t*             result;
    };

    virtual AsyncFile::Reader* getReader();

private:
    int         fd;
};

    PosixAsyncFile*
PosixAsyncFile::open(
    const char* filename,
    bool write)
{
    int fd = ::open(filename, write ? O_CREAT | O_RDWR | O_TRUNC : O_RDONLY, write ? S_IRWXU | S_IRGRP : 0);
    if (fd < 0) {
        WriteErrorMessage("Unable to create SAM file '%s', %d\n",filename,errno);
        return NULL;
    }
    return new PosixAsyncFile(fd);
}

PosixAsyncFile::PosixAsyncFile(
    int i_fd)
    : fd(i_fd)
{
}

    bool
PosixAsyncFile::close()
{
    return ::close(fd) == 0;
}

    AsyncFile::Writer*
PosixAsyncFile::getWriter()
{
    return new Writer(this);
}

PosixAsyncFile::Writer::Writer(PosixAsyncFile* i_file)
    : file(i_file), writing(false)
{
    memset(&aiocb, 0, sizeof(aiocb));
    if (! CreateSingleWaiterObject(&ready)) {
        WriteErrorMessage("PosixAsyncFile: cannot create waiter\n");
        soft_exit(1);
    }
}

    bool
PosixAsyncFile::Writer::close()
{
    waitForCompletion();
    DestroySingleWaiterObject(&ready);
    return true;
}

    void
sigev_ready(
    union sigval val)
{
    SignalSingleWaiterObject((SingleWaiterObject*) val.sival_ptr);
}

    void
aio_setup(
    struct aiocb* control,
    SingleWaiterObject* ready,
    int fd,
    void* buffer,
    size_t length,
    size_t offset)
{
    control->aio_fildes = fd;
    control->aio_buf = buffer;
    control->aio_nbytes = length;
    control->aio_offset = offset;
    control->aio_sigevent.sigev_notify = SIGEV_THREAD;
    control->aio_sigevent.sigev_value.sival_ptr = ready;
    control->aio_sigevent.sigev_notify_function = sigev_ready;
}


    bool
PosixAsyncFile::Writer::beginWrite(
    void* buffer,
    size_t length,
    size_t offset,
    size_t *bytesWritten)
{
    if (! waitForCompletion()) {
        return false;
    }
    aio_setup(&aiocb, &ready, file->fd, buffer, length, offset);
    result = bytesWritten;
    if (aio_write(&aiocb) < 0) {
        warn("PosixAsyncFile aio_write failed");
        return false;
    }
    writing = true;
    return true;
}

    bool
PosixAsyncFile::Writer::waitForCompletion()
{
    if (writing) {
        WaitForSingleWaiterObject(&ready);
	ResetSingleWaiterObject(&ready);
        writing = false;
        ssize_t ret = aio_return(&aiocb);
        if (ret < 0 && errno != 0) {
            warn("PosixAsyncFile Writer aio_return failed");
            return false;
        }
        if (result != NULL) {
            *result = max((ssize_t)0, ret);
        }
    }
    return true;
}

    AsyncFile::Reader*
PosixAsyncFile::getReader()
{
    return new Reader(this);
}

PosixAsyncFile::Reader::Reader(
    PosixAsyncFile* i_file)
    : file(i_file), reading(false)
{
    memset(&aiocb, 0, sizeof(aiocb));
    if (! CreateSingleWaiterObject(&ready)) {
        WriteErrorMessage("PosixAsyncFile cannot create waiter\n");
        soft_exit(1);
    }
}

    bool
PosixAsyncFile::Reader::close()
{
    DestroySingleWaiterObject(&ready);
    return true;
}

    bool
PosixAsyncFile::Reader::beginRead(
    void* buffer,
    size_t length,
    size_t offset,
    size_t* bytesRead)
{
    if (! waitForCompletion()) {
        return false;
    }
    aio_setup(&aiocb, &ready, file->fd, buffer, length, offset);
    result = bytesRead;
    if (aio_read(&aiocb) < 0) {
        warn("PosixAsyncFile Reader aio_read failed");
        return false;
    }
    reading = true;
    return true;
}

    bool
PosixAsyncFile::Reader::waitForCompletion()
{
    if (reading) {
        WaitForSingleWaiterObject(&ready);
	ResetSingleWaiterObject(&ready);
        reading = false;
        ssize_t ret = aio_return(&aiocb);
        if (ret < 0 && errno != 0) {
            warn("PosixAsyncFile Reader aio_return");
            return false;
        }
        if (result != NULL) {
            *result = max((ssize_t)0, ret);
        }
    }
    return true;
}

#else

// todo: make this actually async!

class OsxAsyncFile : public AsyncFile
{
public:
    static OsxAsyncFile* open(const char* filename, bool write);

    OsxAsyncFile(int i_fd);

    virtual bool close();

    class Writer : public AsyncFile::Writer
    {
    public:
        Writer(OsxAsyncFile* i_file);

        virtual bool close();

        virtual bool beginWrite(void* buffer, size_t length, size_t offset, size_t *bytesWritten);

        virtual bool waitForCompletion();
    
    private:
        OsxAsyncFile*       file;
        bool                writing;
        SingleWaiterObject  ready;
        struct aiocb        aiocb;
        size_t*             result;
    };

    virtual AsyncFile::Writer* getWriter();
    
    class Reader : public AsyncFile::Reader
    {
    public:
        Reader(OsxAsyncFile* i_file);

        virtual bool close();

        virtual bool beginRead(void* buffer, size_t length, size_t offset, size_t *bytesRead);

        virtual bool waitForCompletion();
    
    private:
        OsxAsyncFile*       file;
        bool                reading;
        size_t*             result;
    };

    virtual AsyncFile::Reader* getReader();

private:
    int         fd;
};

    OsxAsyncFile*
OsxAsyncFile::open(
    const char* filename,
    bool write)
{
    int fd = ::open(filename, write ? O_CREAT | O_RDWR | O_TRUNC : O_RDONLY, write ? S_IRWXU | S_IRGRP : 0);
    if (fd < 0) {
        WriteErrorMessage("Unable to create SAM file '%s', %d\n",filename,errno);
        return NULL;
    }
    return new OsxAsyncFile(fd);
}

OsxAsyncFile::OsxAsyncFile(
    int i_fd)
    : fd(i_fd)
{
}

    bool
OsxAsyncFile::close()
{
    return ::close(fd) == 0;
}

    AsyncFile::Writer*
OsxAsyncFile::getWriter()
{
    return new Writer(this);
}

OsxAsyncFile::Writer::Writer(OsxAsyncFile* i_file)
    : file(i_file), writing(false)
{
}

    bool
OsxAsyncFile::Writer::close()
{
    return true;
}

    bool
OsxAsyncFile::Writer::beginWrite(
    void* buffer,
    size_t length,
    size_t offset,
    size_t *bytesWritten)
{
    size_t m = ::lseek(file->fd, offset, SEEK_SET);
    if (m == -1) {
        return false;
    }
    size_t n = ::write(file->fd, buffer, length);
    if (bytesWritten) {
      *bytesWritten = n;
    }
    return n != -1;
}

    bool
OsxAsyncFile::Writer::waitForCompletion()
{
    return true;
}

    AsyncFile::Reader*
OsxAsyncFile::getReader()
{
    return new Reader(this);
}

OsxAsyncFile::Reader::Reader(
    OsxAsyncFile* i_file)
    : file(i_file), reading(false)
{
}

    bool
OsxAsyncFile::Reader::close()
{
    return true;
}

    bool
OsxAsyncFile::Reader::beginRead(
    void* buffer,
    size_t length,
    size_t offset,
    size_t* bytesRead)
{
    size_t m = ::lseek(file->fd, offset, SEEK_SET);
    if (m == -1) {
        return false;
    }
    size_t n = ::read(file->fd, buffer, length);
    if (bytesRead) {
        *bytesRead = n;
    }
    return n != -1;
}

    bool
OsxAsyncFile::Reader::waitForCompletion()
{
    return true;
}

#endif

int _fseek64bit(FILE *stream, _int64 offset, int origin)
{
#ifdef __APPLE__
    // Apple's file pointers are already 64-bit so just use fseeko.
    fseeko(stream, offset, origin);
#else
    return fseeko64(stream, offset, origin);
#endif
}

FileMapper::FileMapper()
{
    fd = -1;
    initialized = false;
    mapCount = 0;
    pagesize = getpagesize();
}

bool
FileMapper::init(const char *i_fileName)
{
    if (initialized) {
        if (strcmp(fileName, i_fileName)) {
            WriteErrorMessage("FileMapper already initialized with %s, cannot init with %s\n", fileName, i_fileName);
            return false;
        }
        return true;
    }
    fileName = i_fileName;
    fd = open(fileName, O_RDONLY);
    if (fd == -1) {
	    WriteErrorMessage("Failed to open %s\n", fileName);
	    return false;
    }

    struct stat sb;
    int r = fstat(fd, &sb);
    if (r == -1) {
	    WriteErrorMessage("Failed to stat %s\n", fileName);
	    return false;
    }
    fileSize = sb.st_size;

    initialized = true;

    return true;
}

char *
FileMapper::createMapping(size_t offset, size_t amountToMap, void** o_token)
{
    size_t beginRounding = offset % pagesize;

    size_t mapRequestSize = beginRounding + amountToMap;
    //_ASSERT(mapRequestSize % pagesize == 0);
    if (mapRequestSize + offset >= fileSize) {
        mapRequestSize = 0; // Says to just map the whole thing.
    } 

    char* mappedBase = (char *) mmap(NULL, amountToMap + beginRounding, PROT_READ, MAP_SHARED, fd, offset - beginRounding);
    if (mappedBase == MAP_FAILED) {
	    WriteErrorMessage("mmap failed.\n");
	    return NULL;
    }

    int r = madvise(mappedBase, min((size_t) madviseSize, amountToMap + beginRounding), MADV_WILLNEED | MADV_SEQUENTIAL);
    _ASSERT(r == 0);
    lastPosMadvised = 0;

    InterlockedIncrementAndReturnNewValue(&mapCount);
    *o_token = new UnmapToken(mappedBase, amountToMap);
    return mappedBase + beginRounding;
}

void
FileMapper::unmap(void* i_token)
{
    _ASSERT(mapCount > 0);
    if (mapCount > 0) {
        int n = InterlockedDecrementAndReturnNewValue(&mapCount);
        _ASSERT(n >= 0);
        UnmapToken* token = (UnmapToken*) i_token;
        munmap(token->first, token->second);
        delete token;
    }
}

FileMapper::~FileMapper()
{
    _ASSERT(mapCount == 0);
    close(fd);
}

#if 0
void
FileMapper::prefetch(size_t currentRead)
{
    if (currentRead > lastPosMadvised + madviseSize / 2) {
        _uint64 offset = lastPosMadvised + madviseSize;
        _uint64 len = (offset > amountMapped ? 0 : min(amountMapped - offset, (_uint64) madviseSize));
        if (len > 0) {
            // Start reading new range
            int r = madvise(mappedBase + offset, len, MADV_WILLNEED);
            _ASSERT(r == 0);
        }
        if (lastPosMadvised > 0) {
          // Unload the range we had before our current one
          int r = madvise(mappedBase + lastPosMadvised - madviseSize, madviseSize, MADV_DONTNEED);
          _ASSERT(r == 0);
        }
        lastPosMadvised = offset;
    }
}
#endif

void PreventMachineHibernationWhileThisThreadIsAlive()
{
	// Only implemented for Windows
}

void SetToLowSchedulingPriority()
{
    // Only implemented for Windows (the Linux version is per-thread, and I'm too lazy to do it now).
    WriteErrorMessage("The Linux code for running at low priority is not implemented, so SNAP will run at normal priority\n");
}

//
// Linux named pipes are unidirectional, so we need two of them.
//
struct NamedPipe {
    bool    serverSide;
    char *  pipeName;
    FILE *  input;
    FILE *  output;
};

bool createPipe(const char *fullyQualifiedPipeName)
{
	if (mkfifo(fullyQualifiedPipeName, S_IRUSR | S_IWUSR)) {
		if (errno != EEXIST) {
			if (errno == ENOENT) {
				WriteErrorMessage("OpenNamedPipe: unable to create named pipe at path '%s' because a directory in the path doesn't exist.  Please create it or use a different pipe name\n", fullyQualifiedPipeName);
				return false;
			}

			if (errno == EACCES) {
				WriteErrorMessage("OpenNamedPipe: unable to create named pipe at path '%s' because you do not have sufficient permissions.\n", fullyQualifiedPipeName);
				return false;
			}

			if (errno == ENOTDIR) {
				WriteErrorMessage("OpenNamedPipe: a component of your pipe path isn't a directory.  '%s'\n", fullyQualifiedPipeName);
				return false;
			}

			WriteErrorMessage("OpenNamedPipe: unexpectedly failed to create named pipe '%s', errno %d\n", fullyQualifiedPipeName, errno);
			return false;
		}
	}
	return true;
}

FILE *connectPipe(char *fullyQualifiedPipeName, bool forInput) 
{
	FILE *pipeFile = fopen(fullyQualifiedPipeName, forInput ? "r" : "w");
	if (NULL == pipeFile) {
		WriteErrorMessage("OpenNamedPipe: unable to open pipe file '%s', errno %d\n", fullyQualifiedPipeName, errno);
	}
	return pipeFile;
}

char *createFullyQualifiedPipeName(const char *pipeName, bool serverSide, bool forInput)
{
    char *fullyQualifiedPipeName;
	const char *defaultPipeDirectory = "/tmp/";
	const char *pipeDirectory;
	const char *toServer = "-toServer";
	const char *toClient = "-toClient";

	if (pipeName[0] != '/') {
		pipeDirectory = defaultPipeDirectory;
	} else {
		pipeDirectory = "";
	}

	fullyQualifiedPipeName = new char[strlen(pipeDirectory) + strlen(pipeName) + __max(strlen(toServer), strlen(toClient)) + 1];	// +1 for trailing null

	sprintf(fullyQualifiedPipeName, "%s%s%s", pipeDirectory, pipeName, (serverSide == forInput) ? toServer : toClient);

	return fullyQualifiedPipeName;
}

bool connectNamedPipes(NamedPipe *pipe)
{
    char *inputPipeName = createFullyQualifiedPipeName(pipe->pipeName, pipe->serverSide, true);
    char *outputPipeName = createFullyQualifiedPipeName(pipe->pipeName, pipe->serverSide, false);

    if (pipe->input != NULL) {
        fclose(pipe->input);
        pipe->input = NULL;
    }

    if (pipe->output != NULL) {
        fclose(pipe->output);
        pipe->output = NULL;
    }

    //
    // Connecting pipes is synchronous, so the server and client need to connect
    // in opposite order.
    //
    if (pipe->serverSide) {
        signal(SIGPIPE, SIG_IGN);// If the client hits ^C, we'll get this.  Ignore it, let the fwrite fail, and continue
        pipe->input = connectPipe(inputPipeName,  true);
        pipe->output = connectPipe(outputPipeName, false);
	    //
	    // Release any exclusive lock on the toServer pipe that may have been left by a now-dead client
	    //
	    struct flock lock;
	    lock.l_type = F_UNLCK;
	    lock.l_whence = SEEK_SET;
	    lock.l_start = 0;
	    lock.l_len = 1;
	    lock.l_pid = 0;

	    if (fcntl(fileno(pipe->output), F_SETLKW, &lock) < 0) {
	        fprintf(stderr,"Unable to clear named pipe lock, errno %d\n", errno);
	        delete pipe;
	        return NULL;
	    }
    } else {
        pipe->output = connectPipe(outputPipeName, false);
        pipe->input = connectPipe(inputPipeName,  true);
    }

    delete [] inputPipeName;
    delete [] outputPipeName;

    return pipe->input != NULL && pipe->output != NULL;
}


NamedPipe *OpenNamedPipe(const char *pipeName, bool serverSide)
{
    char *inputPipeName = createFullyQualifiedPipeName(pipeName, serverSide, true);
    char *outputPipeName = createFullyQualifiedPipeName(pipeName, serverSide, false);

	NamedPipe *pipe = new NamedPipe;
	pipe->pipeName = new char[strlen(pipeName) + 1];
	strcpy(pipe->pipeName, pipeName);
	pipe->input = pipe->output = NULL;
	pipe->serverSide = serverSide;

	if (serverSide) {
	    if (!createPipe(inputPipeName)) {
	        delete pipe;
	        return NULL;
	    }
	    if (!createPipe(outputPipeName)) {
	        delete pipe;
	        return NULL;
	    }
	}


	delete [] inputPipeName;
	delete [] outputPipeName;

	if (!connectNamedPipes(pipe)) {
	    delete pipe;
	    return NULL;
	}

	if (!serverSide) {
	    //
	    // Take an exclusive lock on the toServer pipe so that only one client is sending at a time.
	    //
	    struct flock lock;
	    lock.l_type = F_WRLCK;
	    lock.l_whence = SEEK_SET;
	    lock.l_start = 0;
	    lock.l_len = 1;
	    lock.l_pid = 0;

	    if (fcntl(fileno(pipe->output), F_SETLKW, &lock) < 0) {
	        fprintf(stderr,"OpenNamedPipe: F_SETLKW failed, errno %d\n", errno);
	        delete pipe;
	        return NULL;
	    }
	}

	return pipe;
}


//
// Our Linux version of named pipe IO sends strings with 4 byte byte counts first.
//
bool ReadFromNamedPipe(NamedPipe *pipe, char *outputBuffer, size_t outputBufferSize)
{
    unsigned int size;

    for (;;) {
        if (1 != fread(&size, sizeof(size), 1, pipe->input)) {
            if (!pipe->serverSide) {
                return false;
            }
            if (!connectNamedPipes(pipe)) {
                return false;
            }
            continue;
        }

        if (size >= outputBufferSize) {
            WriteErrorMessage("Trying to read too big a chunk from named pipe, %d >= %lld\n", size, outputBufferSize);
            return false;
        }

        if (1 != fread(outputBuffer, size, 1, pipe->input)) {
            if (!pipe->serverSide) {
                return false;
            }
            if (!connectNamedPipes(pipe)) {
                return false;
            }
            continue;
        }

        break;
    }

    outputBuffer[size] = '\0';
    return true;
}

bool WriteToNamedPipe(NamedPipe *pipe, const char *stringToWrite)
{
    unsigned int size = (unsigned int)strlen(stringToWrite);

    if (1 != fwrite(&size, sizeof(size), 1, pipe->output)) {
        return false;
    }

    if (1 != fwrite(stringToWrite, size, 1, pipe->output)) {
        return false;
    }

    fflush(pipe->output);

    return true;
}

void CloseNamedPipe(NamedPipe *pipe)
{
	fclose(pipe->input);
	fclose(pipe->output);
	delete pipe;
}

const char *DEFAULT_NAMED_PIPE_NAME = "SNAP";

#endif  // _MSC_VER

AsyncFile* AsyncFile::open(const char* filename, bool write)
{
    if (!strcmp("-", filename) && write) {
        return StdoutAsyncFile::open("-", true);
    }
#ifdef _MSC_VER
    return WindowsAsyncFile::open(filename, write);
#else
#ifdef __linux__
    return PosixAsyncFile::open(filename, write);
#else
    return OsxAsyncFile::open(filename, write);
#endif
#endif
}

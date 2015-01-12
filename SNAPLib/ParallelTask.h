/*++

Module Name:

    ParallelTask.h

Abstract:

    Simple parallel task manager

Authors:

    Ravi Pandya, May 2012

Environment:

    User mode service.

Revision History:

--*/

#pragma once
#include "stdafx.h"
#include "Compat.h"
#include "exit.h"
#include "Error.h"

/*++
    Simple class to handle parallelized algorithms.
    TContext should extend TContextBase, and provide the following methods:
        void initializeThread()
            Called once on main thread after TContext has been assigned from common,
            and threadNum set.
        void runThread()
            Called to run the thread's work until termination.
            May use something like RangeSplitter to get work.
        void finishThread(TContext* common)
            Called once on main thread after all threads have finished,
            to write results back to common.
--*/
    template <class TContext>
class ParallelTask
{
public:

    inline TContext* getCommonContext() { return common; }

    // i_common should have totalThreads & bindToProcessors set
    ParallelTask(TContext* i_common);

    // run all threads until completion, gather results in common
    void run();

    // run all tasks on a separate thread
    void fork();

private:

    // initial & final context
    TContext*   common;

    // array of per-thread contexts
    TContext*   contexts;

    static void threadWorker(void* threadContext);

    static void forkWorker(void* threadContext);
};

/*++
    Base for type parameter to parallel task
--*/
struct TaskContextBase
{
    // should be set before passing to ParallelTask constructor
    int                 totalThreads;
    bool                bindToProcessors;

    // time taken to run in millis
    _int64              time;

    // for internal use:
    int                 threadNum;        // current thread number, 0...totalThreads-1
    SingleWaiterObject *doneWaiter;       // Gets notified when the last thread ends.
    volatile int        runningThreads;
    volatile int       *pRunningThreads;
#ifdef  _MSC_VER
    volatile int       *nThreadsAllocatingMemory;
    EventObject        *memoryAllocationCompleteBarrier;
    bool                useTimingBarrier;
#endif  // _MSC_VER
};



    template <class TContext>
ParallelTask<TContext>::ParallelTask(
    TContext* i_common)
    : common(i_common), contexts(NULL)
{
    _ASSERT(i_common->totalThreads > 0);
}

    template <class TContext>
    void
ParallelTask<TContext>::run()
{
    _int64 start = timeInMillis();
    SingleWaiterObject doneWaiter;
    if (!CreateSingleWaiterObject(&doneWaiter)) {
        WriteErrorMessage( "Failed to create single waiter object for thread completion.\n");
        soft_exit(1);
    }

#ifdef  _MSC_VER
    int nThreadsAllocatingMemory = common->totalThreads;
    EventObject memoryAllocationCompleteBarrier;
    CreateEventObject(&memoryAllocationCompleteBarrier);

    common->nThreadsAllocatingMemory = &nThreadsAllocatingMemory;
    common->memoryAllocationCompleteBarrier = &memoryAllocationCompleteBarrier;
#endif  // _MSC_VER
    common->doneWaiter = &doneWaiter;
    common->runningThreads = common->totalThreads;
    common->pRunningThreads = &common->runningThreads;

    contexts = new TContext[common->totalThreads];
    for (int i = 0; i < common->totalThreads; i++) {
        contexts[i] = *common;
        contexts[i].threadNum = i;
        contexts[i].initializeThread();

        if (!StartNewThread(ParallelTask<TContext>::threadWorker, &contexts[i])) {
            WriteErrorMessage( "Unable to start worker thread.\n");
            soft_exit(1);
        }
    }

#ifdef  _MSC_VER
    if (common->useTimingBarrier) {
        WaitForEvent(&memoryAllocationCompleteBarrier);
        WriteStatusMessage("Cleared timing barrier.\n");
        start = timeInMillis();
    }
#endif  // _MSC_VER

    if (!WaitForSingleWaiterObject(&doneWaiter)) {
        WriteErrorMessage( "Waiting for all threads to finish failed\n");
        soft_exit(1);
    }
    DestroySingleWaiterObject(&doneWaiter);
#ifdef  _MSC_VER
    DestroyEventObject(&memoryAllocationCompleteBarrier);
#endif  // _MSC_VER

    for (int i = 0; i < common->totalThreads; i++) {
        contexts[i].finishThread(common);
    }

    common->time = timeInMillis() - start;
}

    template <class TContext>
    void
ParallelTask<TContext>::fork()
{
    if (!StartNewThread(ParallelTask<TContext>::forkWorker, this)) {
        WriteErrorMessage( "Unable to fork task thread.\n");
        soft_exit(1);
    }
}

    template <class TContext>
    void
ParallelTask<TContext>::forkWorker(
    void* forkArg)
{
    ((ParallelTask<TContext>*) forkArg)->run();
}

    template <class TContext>
    void
ParallelTask<TContext>::threadWorker(
    void* threadArg)
{
    TContext* context = (TContext*) threadArg;
    if (context->bindToProcessors) {
        BindThreadToProcessor(context->threadNum);
    }

    context->runThread();

    // Decrement the running thread count and wake up the waiter if it hits 0.
    if (0 == InterlockedDecrementAndReturnNewValue(context->pRunningThreads)) {
        SignalSingleWaiterObject(context->doneWaiter);
	}
}

struct WorkerContext;
class ParallelWorker;
class ParallelWorkerManager;

// coroutined parallel workers
// does code inline if numThreads = 0
// can either callback when done, or synchronously wait for all to complete
class ParallelCoworker
{
public:

    typedef void (*Callback)(void*);
    ParallelCoworker(int i_numThreads, bool i_bindToProcessors, ParallelWorkerManager* supplier, Callback callback = NULL, void* parameter = NULL);

    ~ParallelCoworker();

    // start forked thread running
    void start();

    // do one unit of work, asynchronously if callback, else synchronously
    void step();

    // stop everything, waits until all threads exit
    void stop();

    ParallelWorkerManager* getManager() { return manager; }

private:
    EventObject *workReady; // One per worker thread
    EventObject *workDone;  // One per worker thread
    ParallelWorker** workers; // one per worker thread
    ParallelWorkerManager* manager;
    volatile bool stopped;
    const int numThreads;
    const bool bindToProcessors;
    Callback callback;
    void* parameter;
    WorkerContext* context;
    ParallelTask<WorkerContext>* task;
    SingleWaiterObject finished;

    friend struct WorkerContext;
};

// abstract classes for specifying the actual work

// creates new per-thread workers
class ParallelWorker;

class ParallelWorkerManager
{
public:

    // todo: using void* context pointers to avoid pain of templates but should really be made typesafe
    virtual void initialize(void* context) {}

    virtual ParallelWorker* createWorker() = 0;

    virtual void beginStep() {}

    virtual void finishStep() {}

    void configure(ParallelWorker* worker, int threadNum, int totalThreads); // special case
};

// per-thread worker
class ParallelWorker
{
public:
    ParallelWorker() {}
    virtual ~ParallelWorker() {}

    virtual void initialize() {}
    virtual void step() = 0;

protected:
    ParallelWorkerManager* getManager() { return manager; }
    int getThreadNum() { return threadNum; }
    int getNumThreads() { return numThreads; }

private:
    friend class ParallelCoworker;
    friend class ParallelWorkerManager;
    void configure(ParallelWorkerManager* i_manager, int i_threadNum, int i_numThreads)
    { manager = i_manager; threadNum = i_threadNum; numThreads = i_numThreads; }

    ParallelWorkerManager* manager;
    int threadNum;
    int numThreads;
};

struct WorkerContext : public TaskContextBase
{
    ParallelCoworker* shared;

    void initializeThread();
    void runThread();
    void finishThread(WorkerContext* common);
};


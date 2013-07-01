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
        fprintf(stderr, "Failed to create single waiter object for thread completion.\n");
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
            fprintf(stderr, "Unable to start worker thread.\n");
            soft_exit(1);
        }
    }

#ifdef  _MSC_VER
    if (common->useTimingBarrier) {
        WaitForEvent(&memoryAllocationCompleteBarrier);
        printf("Cleared timing barrier.\n");
        start = timeInMillis();
    }
#endif  // _MSC_VER

    if (!WaitForSingleWaiterObject(&doneWaiter)) {
        fprintf(stderr, "Waiting for all threads to finish failed\n");
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
        fprintf(stderr, "Unable to fork task thread.\n");
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

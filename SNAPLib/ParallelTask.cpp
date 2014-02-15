/*++

Module Name:

    ParallelTask.cpp

Abstract:

    Parallel task management

Environment:

    User mode service.

--*/

#include "stdafx.h"
#include "ParallelTask.h"
#include "Error.h"

using std::max;

ParallelCoworker::ParallelCoworker(int i_numThreads, bool i_bindToProcessors, ParallelWorkerManager* i_manager, Callback i_callback, void* i_parameter)
    : stopped(false), numThreads(i_numThreads), bindToProcessors(i_bindToProcessors), manager(i_manager), callback(i_callback), parameter(i_parameter)
{
    workReady = new EventObject[numThreads];
    workDone = new EventObject[numThreads];
    workers = new ParallelWorker*[max(numThreads, 1)];
    for (int i = 0; i < numThreads; i++) {
        CreateEventObject(&workReady[i]);
        PreventEventWaitersFromProceeding(&workReady[i]);
        CreateEventObject(&workDone[i]);
        PreventEventWaitersFromProceeding(&workDone[i]);
        workers[i] = manager->createWorker();
        workers[i]->configure(manager, i, numThreads);
    }
    CreateSingleWaiterObject(&finished);
}

ParallelCoworker::~ParallelCoworker()
{
    for (int i = 0; i < numThreads; i++) {
        DestroyEventObject(&workReady[i]);
        DestroyEventObject(&workDone[i]);
        delete workers[i];
    }
    delete [] workReady;
    delete [] workDone;
    delete [] workers;
    delete task;
    delete context;
}

void ParallelCoworker::start()
{
    context = new WorkerContext();
    context->shared = this;
    context->totalThreads = numThreads;
    context->bindToProcessors = bindToProcessors;
#ifdef _MSC_VER
    context->useTimingBarrier = false;
#endif
    task = new ParallelTask<WorkerContext>(context);
    task->fork();
}

void ParallelCoworker::step()
{
    manager->beginStep();
    for (int i = 0; i < numThreads; i++) {
        PreventEventWaitersFromProceeding(&workDone[i]);
        AllowEventWaitersToProceed(&workReady[i]);
    }
    // if async, thread 0 will callback when all workers finish
    // if sync, wait for all workers to finish
    if (callback == NULL) {
        for (int i = 0; i < numThreads; i++) {
            WaitForEvent(&workDone[i]);
        }
        manager->finishStep();
    }
}

void ParallelCoworker::stop()
{
    stopped = true;
    for (int i = 0; i < numThreads; i++) {
        AllowEventWaitersToProceed(&workReady[i]);
    }
    if (!WaitForSingleWaiterObject(&finished)) {
        WriteErrorMessage("Waiting for all threads to finish failed\n");
        soft_exit(1);
    }
}

    void
WorkerContext::initializeThread()
{
    shared->workers[threadNum]->initialize();
}

    void
WorkerContext::runThread()
{
    while (true) {
        //fprintf(stderr, "worker task thread %d waiting to begin\n", GetCurrentThreadId());
        WaitForEvent(&shared->workReady[threadNum]);
        PreventEventWaitersFromProceeding(&shared->workReady[threadNum]);
        if (shared->stopped) {
            return;
        }
        //fprintf(stderr, "worker task thread %d begin\n", GetCurrentThreadId());
        _int64 start = timeInMillis();
        shared->workers[threadNum]->step();
        //fprintf(stderr, "zip task thread %d done %lld ms\n", GetCurrentThreadId(), timeInMillis() - start);
        AllowEventWaitersToProceed(&shared->workDone[threadNum]);
        
        // if async, thread 0 will wait for everyone to finish and invoke callback
        if (threadNum == 0 && shared->callback != NULL) {
            for (int i = 1; i < shared->numThreads; i++) {
                WaitForEvent(&shared->workDone[i]);
            }
            shared->manager->finishStep();
            shared->callback(shared->parameter);
        }
    }
}

    void
WorkerContext::finishThread(WorkerContext* common)
{
    if (threadNum == totalThreads - 1) {
        SignalSingleWaiterObject(&shared->finished);
    }
}

    void
ParallelWorkerManager::configure(
    ParallelWorker* worker,
    int threadNum,
    int totalThreads)
{
    worker->configure(this, threadNum, totalThreads);
}

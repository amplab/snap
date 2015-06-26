#include "stdafx.h"
#include "TestLib.h"
#include "LandauVishkin.h"

struct EventTest;

struct TestContext {
  EventTest* parent;
  int index;
  SingleWaiterObject singleWaiterObject;

  void init(EventTest* parent, int i);
};

struct EventTest {
  EventObject event;
  TestContext* contexts;
  bool bind;
  volatile int started;
  volatile int proceeded;

  void testManyWaiters(int threads, bool bind);
  void testSingleWaiters(int threads, bool bind);
};

void TestContext::init(EventTest* i_parent, int i) {
    this->parent = i_parent;
    index = i;
    CreateSingleWaiterObject(&singleWaiterObject);
}

void waitEqual(int value, volatile int* variable, const char* message)
{
  for (int i = 0; i < 1000; i++) {
    if (*variable == value) {
      return;
    }
    SleepForMillis(5);
  }
  ASSERT_EQ_M(value, *variable, message);
}

void testManyWaitersMain(void* c)
{
  TestContext* context = (TestContext*) c;
  if (context->parent->bind) {
    BindThreadToProcessor(context->index%GetNumberOfProcessors());
  }
  //printf("start thread %d%s\n", context->index, context->parent->bind ? " bind" : "");
  InterlockedIncrementAndReturnNewValue(&context->parent->started);
  WaitForEvent(&context->parent->event);
  //printf("proceed thread %d\n", context->index);
  InterlockedIncrementAndReturnNewValue(&context->parent->proceeded);
}

void EventTest::testManyWaiters(int threads, bool i_bind)
{
  //printf("testing many %d threads%s\n", threads, i_bind ? " (bind)" : "");
  contexts = new TestContext[threads];
  CreateEventObject(&event);
  PreventEventWaitersFromProceeding(&event);
  started = proceeded = 0;
  bind = i_bind;
  for (int i = 0; i < threads; i++) {
    contexts[i].init(this, i);
    StartNewThread(testManyWaitersMain, &contexts[i]);
  }
  char buf[100];
  sprintf(buf, "many started %d%s\n", threads, i_bind ? " bind" : "");
  waitEqual(threads, &started, buf);
  AllowEventWaitersToProceed(&event);
  sprintf(buf, "many proceeded %d%s\n", threads, i_bind ? " bind" : "");
  waitEqual(threads, &proceeded, buf);
  DestroyEventObject(&event);
}

TEST_F(EventTest, "many waiters") {
  for (int i = 0; i < 100; i++) {
    for (int threads = 1; threads <= 64; threads *= 2) {
      for (int bind = 0; bind < (threads <= 16 ? 2 : 1); bind++) {
	testManyWaiters(threads, bind);
      }
    }
  }
}

void testSingleWaitersMain(void* c)
{
  TestContext* context = (TestContext*) c;
  if (context->parent->bind) {
    BindThreadToProcessor(context->index%GetNumberOfProcessors());
  }
  //printf("start thread %d%s\n", context->index, context->parent->bind ? " bind" : "");
  InterlockedIncrementAndReturnNewValue(&context->parent->started);
  WaitForSingleWaiterObject(&context->singleWaiterObject);
  //printf("proceed thread %d\n", context->index);
  InterlockedIncrementAndReturnNewValue(&context->parent->proceeded);
  ResetSingleWaiterObject(&context->singleWaiterObject);
}

void EventTest::testSingleWaiters(int threads, bool i_bind)
{
  //printf("testing single %d threads%s\n", threads, i_bind ? " (bind)" : "");
  contexts = new TestContext[threads];
  started = proceeded = 0;
  bind = i_bind;
  for (int i = 0; i < threads; i++) {
    contexts[i].init(this, i);
    StartNewThread(testSingleWaitersMain, &contexts[i]);
  }
  //SleepForMillis(50);
  char buf[100];
  sprintf(buf, "single started %d%s\n", threads, i_bind ? " bind" : "");
  waitEqual(threads, &started, buf);
  for (int i = 0; i < threads; i++) {
    SignalSingleWaiterObject(&contexts[i].singleWaiterObject);
    sprintf(buf, "single proceeded %d of %d%s\n", i, threads, i_bind ? " bind" : "");
    waitEqual(i + 1, &proceeded, buf);
    //SleepForMillis(10); // allow more threads to release
    sprintf(buf, "single after proceeded %d of %d%s\n", i, threads, i_bind ? " bind" : "");
    ASSERT_EQ_M(i + 1, proceeded, buf);
  }

  for (int i = 0; i < threads; i++) {
      DestroySingleWaiterObject(&contexts[i].singleWaiterObject);
  }
}

TEST_F(EventTest, "single waiters") {
  for (int i = 0; i < 50; i++) {
    for (int threads = 1; threads <= 64; threads *= 2) {
      for (int bind = 0; bind < (threads <= 16 ? 2 : 1); bind++) {
	    testSingleWaiters(threads, bind); 
      }
    }
  }
}

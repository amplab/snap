#pragma once

#include "Compat.h"
#include "stdafx.h"
#include "BigAlloc.h"
#include "VariableSizeMap.h"

using std::max;
using std::min;

// adapted from http://visualstudiomagazine.com/articles/2012/11/01/priority-queues-with-c.aspx

template <typename P, typename V>
class PriorityQueue
{
private:
    struct Entry
    {
        Entry() : value(), priority(0) {}
        Entry(V v, P p) : value(v), priority(p) {}
        void operator=(const Entry& e) { value = e.value; priority = e.priority; }
        V value;
        P priority;
    };
    typedef VariableSizeVector<Entry> EntryVector;

#if 0
    inline void check()
    { _ASSERT(validate()); }
#else
    inline void check() {}
#endif

public:
    // add an element with specific priority
    void add(V value, P pri)
    {
        queue.push_back(Entry(value, pri));
        _int64 ci = queue.size() - 1;
        while (ci > 0) {
            _int64 pi = (ci - 1) / 2; // parent index
            if (queue[ci].priority >= queue[pi].priority) {
                break; // child >= parent so stop
            }
            Entry tmp = queue[ci]; queue[ci] = queue[pi]; queue[pi] = tmp;
            ci = pi;
        }
        check();
    }

    V pop(P* o_priority = NULL)
    {
        _int64 li = queue.size() - 1; // last index (before removal)
        Entry frontItem = queue[0];   // fetch the front
        queue[0] = queue[li];
        queue.erase(li);

        --li; // last index (after removal)
        int pi = 0; // parent index. start at front of pq
        while (true) {
            int ci = pi * 2 + 1; // left child index of parent
            if (ci > li) {
                break;  // no children so done
            }
            int rc = ci + 1;     // right child
            if (rc <= li && queue[rc].priority < queue[ci].priority) {
                // if there is a rc (ci + 1), and it is smaller than left child, use the rc instead
                ci = rc;
            }
            if (queue[pi].priority <= queue[ci].priority) {
                break; // parent is smaller than (or equal to) smallest child so done
            }
            Entry tmp = queue[pi]; queue[pi] = queue[ci]; queue[ci] = tmp; // swap parent and child
            pi = ci;
        }
        check();
        if (o_priority != NULL) {
            *o_priority = frontItem.priority;
        }
        return frontItem.value;
    }

    V peek(P* o_priority = NULL) const
    {
        if (o_priority != NULL) {
            *o_priority = queue[0].priority;
        }
        return queue[0].value;
    }

    void clear()
    { queue.clear(); }

    _int64 size() const
    { return queue.size(); }

    bool validate()
    {
        // is the heap property true for all queue?
        if (queue.size() == 0) {
            return true;
        }
        int li = queue.size() - 1; // last index
        for (int pi = 0; pi <= li; ++pi) {
            int lci = 2 * pi + 1; // left child index
            int rci = 2 * pi + 2; // right child index

            if (lci <= li && queue[pi].priority > queue[lci].priority) return false; // if lc exists and it's greater than parent then bad.
            if (rci <= li && queue[pi].priority > queue[rci].priority) return false; // check the right child too.
        }
        return true;
    }

private:

    EntryVector queue; // sorted list
};

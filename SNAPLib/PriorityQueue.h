#pragma once

#include "Compat.h"
#include "stdafx.h"
#include "BigAlloc.h"
#include "VariableSizeMap.h"

using std::max;
using std::min;

template <typename P, typename V, int _empty=0>
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
        static bool comparator(const Entry& a, const Entry& b)
        { return a.priority > b.priority; }
    };
    typedef VariableSizeVector<Entry> EntryVector;

public:
    // add an element with specific priority, or set priority if exists
    // todo: array is not the most efficient implementation, but it's easy
    // return true if added, false if exists
    bool put(V value, P pri)
    {
        P* p = priority.tryFind(value);
        if (p == NULL) {
            priority.put(value, pri);
            Entry e(value, pri);
            queue.insert(e, Entry::comparator);
            return true;
        } else if (*p != pri) {
            // find where it is in the queue
            P old = *p;
            Entry e(value, old);
            typename EntryVector::iterator i = std::lower_bound(queue.begin(), queue.end(), e, Entry::comparator);
            _ASSERT(i != queue.end());
            while (i->value != value) {
                i++;
                _ASSERT(i != queue.end());
            }
            Entry e2(value, pri);
            if (pri > old) {
                // move towards front
                Entry* j = std::upper_bound(queue.begin(), i, e2, Entry::comparator);
                if (j < i) {
                    memmove(j + 1, j, (i - j) * sizeof(Entry));
                }
                *j = e2;
            } else {
                // move towards back
                typename EntryVector::iterator j = std::lower_bound(i, queue.end(), e2, Entry::comparator);
                if (j == queue.end()) {
                    _ASSERT(i == j - 1);
                    j = i;
                } else if (i < j) {
                    memmove(i, i + 1, (j - i) * sizeof(Entry));
                }
                *j = e2;
            }
            *p = pri;
        }
#if 1
        _ASSERT(queue.size() == priority.size());
        for (Entry* i = queue.begin(); i != queue.end(); i++) {
            _ASSERT(priority[i->value] == i->priority);
            if (i != queue.begin()) {
                _ASSERT(i->priority <= (i-1)->priority);
            }
        }
#endif
        return false;
    }

    // increment priority of element, or add with priority 1 if not there, return new priority
    int increment(V value)
    {
        P* p = priority.tryFind(value);
        P n = p ? *p + 1 : 1;
        put(value, n);
        return n;
    }

    // check if an element is in the list
    bool contains(V value) const
    { return priority.tryFind(value) != NULL; }

    // remove an element with the highest priority from the queue
    V pop(P* o_priority = NULL)
    {
        Entry result = queue[0];
        queue.remove(queue.begin());
        priority.erase(result.value);
        if (o_priority != NULL) {
            *o_priority = result.priority;
        }
        return result.value;
    }

    V peek(P* o_priority = NULL) const
    {
        if (o_priority != NULL) {
            *o_priority = queue[0].priority;
        }
        return queue[0].value;
    }

    // clear all elements
    void clear()
    {
        priority.clear();
        queue.clear();
    }

    int size() const
    { return queue.size(); }

private:

    VariableSizeMap<V,P,150,MapNumericHash<V>,90,_empty> priority; // value -> priority
    EntryVector queue; // sorted list
};

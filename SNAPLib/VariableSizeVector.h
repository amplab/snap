#pragma once

#include "Util.h"

//
// A variable-size vector that does not perform any memory allocation except to grow.
//
template<typename V, int grow = 150, bool big = false>
class VariableSizeVector
{
private:

    inline static void* allocate(size_t bytes)
    {
#ifdef USE_DEVTEAM_OPTIONS
        if (bytes > (1L << 23) && ! big) {
            WriteErrorMessage("%s: allocate %lld - consider using BigAlloc\n", __FUNCTION__, bytes);
        }
#endif
        return big ? BigAlloc(bytes) : malloc(bytes);
    }

    inline static void deallocate(void* p)
    {
        if (big) { BigDealloc(p); } else { free(p); }
    }

public:
    VariableSizeVector(int i_capacity = 16)
        : entries(NULL), count(0), capacity(i_capacity)
    {}
    
    VariableSizeVector(VariableSizeVector& other)
        : entries(other.entries), count(other.count), capacity(other.capacity)
    {
        other.count = 0;
        other.entries = NULL;
    }

    ~VariableSizeVector()
    {
        if (entries != NULL) {
            deallocate(entries);
            entries = NULL;
            count = 0;
        }
    }
    
private:
    inline void increase()
    {
        if (entries == NULL) {
            reserve(capacity);
        } else if (count == capacity) {
            reserve((int) (((_int64) count * grow) / 100));
        }
    }

public:
    void operator=(VariableSizeVector<V>& other)
    {
        entries = other.entries;
        capacity = other.capacity;
        count = other.count;
        other.entries = NULL;
        other.count = 0;
    }

    void reserve(_int64 newCapacity)
    {
        _ASSERT(newCapacity >= 0);
        if (newCapacity <= capacity && entries != NULL) {
            return;
        }
        V* old = entries;
        capacity = __max(newCapacity, capacity);
        entries = (V*) allocate(capacity * sizeof(V));
        if (old != NULL) {
            memcpy(entries, old, count * sizeof(V));
            deallocate(old);
        }
    }

    inline void clear()
    {
        count = 0;
    }

    inline void clean()
    {
        if (entries != NULL) {
            deallocate(entries);
            entries = NULL;
            count = 0;
        }
    }

    inline _int64 size() const
    {
        return count;
    }

    void truncate(int newCount)
    {
		if (newCount < count) {
			count = newCount;
		}
    }
    
    inline void push_back(V& value)
    {
        if (entries == NULL) {
            reserve(capacity);
        } else if (count == capacity) {
            reserve((int) (((_int64) count * grow) / 100));
        }
        _ASSERT(count < capacity);
        entries[count++] = value;
    }
    
    inline void push_back(const V& value)
    {
        increase();
        _ASSERT(count < capacity);
        entries[count++] = value;
    }

    inline void append(VariableSizeVector<V>* other)
    {
        if (other->count == 0) {
            return;
        }
        reserve(count + other->count);
        // todo: allow for operator assign/copy constructor?
        memcpy(&entries[count], other->entries, other->count * sizeof(V));
        count += other->count;
    }
    
    typedef bool comparator(const V& a, const V& b);

    inline int insertionIndex(const V& value, comparator compare, bool before = false)
    {
        V* p = before ? std::lower_bound(entries, entries + count, value, compare)
            : std::upper_bound(entries, entries + count, value, compare);
        int index = (int) (p - entries);
        _ASSERT(index >= 0 && index <= count);
        return index;
    }
    
    // insert into sorted list, AFTER existing elements with same value
    inline int insert(const V& value, comparator compare, bool before = false)
    {
        int index = insertionIndex( value, compare, before);
        increase(); // todo: could fold memmove into new array copy to save time...
        _ASSERT(count < capacity);
        if (index < count) {
            memmove(entries + (index + 1), entries + index, (count - index) * sizeof(V));
        }
        entries[index] = value;
        count++;
        return index;
    }

    inline bool add(const V& value)
    {
        for (int i = 0; i < count; i++) {
            if (entries[i] == value) {
                return false;
            }
        }
        push_back(value);
        return true;
    }

    inline void erase(_int64 index)
    {
        _ASSERT(index >= 0 && index < count);
        if (index < 0 || index >= count) {
            return;
        }
        if (index < count - 1) {
            memmove(entries + index, entries + index + 1, (count - index - 1) * sizeof(V));
        }
        count--;
    }

    inline void extend(int size)
    {
        if (count < size) {
            reserve(size);
            memset(((V*)entries) + count, 0, sizeof(V) * (size - count));
            count = size;
        }
    }

    inline V& operator[](_int64 index) const
    {
        _ASSERT(index >= 0 && index < count);
        return entries[index];
    }

    typedef V* iterator;

    inline iterator findRange(const V& low, const V& high, comparator compare, iterator* o_end)
    {
        *o_end = entries + insertionIndex(high, compare, true);
        return entries + insertionIndex(low, compare, true);
    }
    
    // unsorted search
    inline iterator search(const V& value)
    {
        for (iterator i = begin(); i != end(); i++) {
            if (*i == value) {
                return i;
            }
        }
        return end();
    }
    iterator begin()
    {
        return entries;
    }

    iterator end()
    {
        return &entries[count];
    }
    
    inline void remove(iterator p)
    {
        _ASSERT(p >= entries && p < entries + count);
        if (p < entries + count - 1) {
            memmove(p, p + 1, (count - (p - entries) - 1) * sizeof(V));
        }
        count--;
    }

private:
    V *entries;
    _int64 capacity;
    _int64 count;
};

typedef VariableSizeVector<unsigned> IdVector;
typedef VariableSizeVector<int> IntVector;
using util::IdPair;
using util::IdIntPair;
typedef VariableSizeVector<IdPair> IdPairVector;
typedef VariableSizeVector<IdIntPair> IdIntPairVector;

#pragma once

//
// A variable-size vector that does not perform any memory allocation except to grow.
//
template<typename V, int grow = 150, bool big = false>
class VariableSizeVector
{
private:

    inline static void* allocate(size_t bytes)
    {
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

    void operator=(VariableSizeVector<V>& other)
    {
        entries = other.entries;
        capacity = other.capacity;
        count = other.count;
        other.entries = NULL;
        other.count = 0;
    }

    void reserve(int newCapacity)
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

    inline int size()
    {
        return count;
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
        if (entries == NULL) {
            reserve(capacity);
        } else if (count == capacity) {
            reserve((int) (((_int64) count * grow) / 100));
        }
        _ASSERT(count < capacity);
        entries[count++] = value;
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

    inline void erase(int index)
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

    inline V& operator[](int index)
    {
        _ASSERT(index >= 0 && index < count);
        return entries[index];
    }

    typedef V* iterator;

    iterator begin()
    {
        return entries;
    }

    iterator end()
    {
        return &entries[count];
    }

private:
    V *entries;
    int capacity;
    int count;
};

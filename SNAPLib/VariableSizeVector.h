#pragma once

//
// A variable-size vector that does not perform any memory allocation except to grow.
//
template<typename V, int grow = 150>
class VariableSizeVector
{
public:
    VariableSizeVector(int i_capacity = 16)
        : entries(NULL), count(0), capacity(i_capacity)
    {
    }

    ~VariableSizeVector()
    {
        if (entries != NULL) {
            delete[] entries;
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

    void reserve(int capacity)
    {
        V* old = entries;
        this->capacity = capacity;
        entries = new V[capacity];
        if (old != NULL) {
            memcpy(entries, old, count * sizeof(V));
            delete[] old;
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

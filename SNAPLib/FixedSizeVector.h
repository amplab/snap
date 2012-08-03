#pragma once

//
// A fixed-size vector that does not perform any memory allocation.
//
template<typename V>
class FixedSizeVector
{
public:
    FixedSizeVector(int capacity_ = 16): entries(NULL), curSize(0) {
        reserve(capacity_);
    }

    // Create a fixed size vector initialized to size copies of an initialValue
    FixedSizeVector(int size, V initialValue): entries(NULL), curSize(0) {
        reserve(size);
        for (int i = 0; i < size; i++) {
            push_back(initialValue);
        }
    }

    ~FixedSizeVector() {
        delete[] entries;
    }

    void reserve(int capacity) {
        if (entries != NULL) {
            if (curSize > 0) {
                fprintf(stderr, "reserve() called on a non-empty FixedSizeVector\n");
                exit(1);
            }
            delete[] entries;
        }
        this->capacity = capacity;
        entries = new V[capacity];
    }

    void clear() {
        curSize = 0;
    }

    int size() {
        return curSize;
    }

    inline void push_back(const V& value) {
        _ASSERT(curSize < capacity);
        entries[curSize++] = value;
    }

    inline V& operator[] (int index) {
        return entries[index];
    }

    typedef V *iterator;

    iterator begin() {
        return entries;
    }

    iterator end() {
        return entries + curSize;
    }

private:
    V *entries;
    int capacity;
    int curSize;
};

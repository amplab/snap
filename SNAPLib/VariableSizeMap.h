#pragma once

#include "Compat.h"
#include "BigAlloc.h"

//
// A hash function for numeric types.
//
template<typename T>
class MapNumericHash
{
public:
    inline _uint64 operator() (T value) const {
        return  ((_uint64)value * 131);
    }
};

template<typename K, typename V>
struct VariableSizeMapEntry
{
    VariableSizeMapEntry() : key(), value() {}
    VariableSizeMapEntry(K k, V v) : key(k), value(v) {}
    K key;
    V value;
};

// 
// A variable-size hash map that allows automatic growth
// and does not perform any memory allocation except when growing.
// Allows multi-threaded put, as long as growth=0 (i.e. fixed-size)
// Shared base class for single- and multi-valued maps
//
using std::max;
using std::min;
template<
    typename K,
    typename V,
    int growth = 150,
    typename Hash = MapNumericHash<K>,
    int fill = 80,
    int _empty = 0,
    int _tombstone = -1,
    bool multi = false,
    bool _big = false>
class VariableSizeMapBase
{
protected:
    VariableSizeMapBase(int i_capacity = 16)
        : entries(NULL), count(0), capacity(i_capacity), occupied(0)
    {
        reserve(max(16,i_capacity));
    }

    VariableSizeMapBase(void** data, unsigned i_capacity)
        : entries((Entry*) (3 + (_int64*) *data)),
        capacity(i_capacity),
        count((int) ((_int64*)*data)[0]),
        limit((int) ((_int64*)*data)[1]),
        occupied((int) ((_int64*)*data)[2])
    {
        *data = ((char*)*data) + (size_t) i_capacity * sizeof(Entry) + 3 * sizeof(_int64);
    }

    inline void grow()
    {
        _ASSERT(growth > 100);
        _int64 larger = ((_int64) capacity * growth) / 100;
        _ASSERT(larger < INT32_MAX);
        reserve((int) larger);
    }

    inline void assign(VariableSizeMapBase<K,V>* other)
    {
        if (entries != NULL) {
            if (_big) {
                BigDealloc(entries);
            } else {
                delete [] entries;
            }
        }
        entries = other->entries;
        capacity = other->capacity;
        count = other->count;
        limit = other->limit;
        hash = other->hash;
        occupied = other->occupied;
        other->entries = NULL;
        other->count = 0;
    }
    
public:
    
    inline int size()
    { return count; }

    inline int getCapacity()
    { return capacity; }
    
    ~VariableSizeMapBase()
    {
        if (entries != NULL) {
            if (_big) {
                BigDealloc(entries);
            } else {
                delete [] entries;
            }
        }
        entries = NULL;
        count = 0;
    }
    
    void reserve(int larger)
    {
        Entry* old = entries;
        int small = capacity;
        capacity = larger;
        if (_big) {
            entries = (Entry*) BigAlloc(larger * sizeof(Entry));
        } else {
            entries = new Entry[larger];
        }
        _ASSERT(entries != NULL);
        clear();
        count = 0;
        // grow before it gets to a certain fraction; always leave 1 slot for empty sentinel
        limit = growth == 0 ? capacity - 1 : min(capacity - 1, (int) (((_int64) capacity * fill) / 100));
        _ASSERT(limit > 0);
        if (old != NULL) {
            for (int i = 0; i < small; i++) {
                K k = old[i].key;
                if (k != _empty && k != _tombstone) {
                    Entry* p = this->scan(k, true);
                    _ASSERT(p != NULL);
                    p->key = k;
                    p->value = old[i].value;
                    count++;
                }
            }
            occupied = count;
            if (_big) {
                BigDealloc(old);
            } else {
                delete [] old;
            }
        }
    }
    
    void clear()
    {
        if (entries != NULL) {
            if (_empty == 0 && sizeof(Entry) < 4 * sizeof(K) && ! _big) {
                // optimize zero case
                memset(entries, 0, capacity * sizeof(Entry));
            } else {
                const K e(_empty);
                for (int i = 0; i < capacity; i++) {
                    entries[i].key = e;
                }
            }
        }
        count = occupied = 0;
    }
    
    typedef VariableSizeMapEntry<K,V> Entry;

    typedef Entry* iterator;

    iterator begin()
    {
        return next(&entries[-1]);
    }

    iterator next(iterator x)
    {
        Entry* final = &entries[capacity];
        if (x < final) {
            do {
                x++;
            } while (x < final && (x->key == _empty || x->key == _tombstone));
        }
        return x;
    }

    iterator end()
    {
        return &entries[capacity];
    }

    iterator find(K key)
    {
        Entry* p = this->scan(key, false);
        return p != NULL ? p : end();
    }

    void writeFile(LargeFileHandle* file)
    {
        _int64 x = (_int64) count;
        WriteLargeFile(file, &x, sizeof(_int64));
        x = (_int64) limit;
        WriteLargeFile(file, &x, sizeof(_int64));
        x = (_int64) occupied;
        WriteLargeFile(file, &x, sizeof(_int64));
        WriteLargeFile(file, entries, sizeof(Entry) * (size_t) capacity);
    }
    
protected:
    
    static const int MaxQuadraticProbes = 3;

    void init(int& pos, int& i, K key)
    {
        _ASSERT(key != _empty && key != _tombstone);
        pos = hash(key) % capacity;
        i = 1;
        if (entries == NULL) {
            reserve(capacity);
        }
    }

    bool advance(int& pos, int& i, K key) const
    {
        if (i >= capacity + MaxQuadraticProbes) {
            pos = capacity;
            return false;
        }
        pos = (pos + (i <= MaxQuadraticProbes ? i : 1)) % capacity;
        i++;
        return true;
    }

    Entry* scan(K key, bool add)
    {
        int pos;
        int i;
        init(pos, i, key);
        if (pos == capacity) {
            return NULL;
        }
        while (true) {
            Entry* p = &entries[pos];
            K k = p->key;
            if (k == key && ! (multi && add)) {
                return p;
            } else if (k == _empty) {
                return add ? p : NULL;
            } else if (add && k == _tombstone && ! multi) {
                return p;
            } else if (! advance(pos, i, key)) {
                return NULL;
            }
        }
    }

    Entry *entries;
    int capacity;
    int count;
    int occupied; // number of non-empty slots (includes tombstones)
    int limit; // current limit (capacity * fill / 100)
    Hash hash;
};

// 
// Single-valued map
//
template< typename K, typename V, int growth = 150, typename Hash = MapNumericHash<K>,
    int fill = 80, int _empty = 0, int _tombstone = -1, bool _big = false >
class VariableSizeMap
    : public VariableSizeMapBase<K,V,growth,Hash,fill,_empty,_tombstone,false,_big>
{
public:
    VariableSizeMap(int i_capacity = 16)
        : VariableSizeMapBase<K,V,growth,Hash,fill,_empty,_tombstone,false,_big>(i_capacity)
    {}

    VariableSizeMap(const VariableSizeMap<K,V>& other)
    {
        this->assign((VariableSizeMapBase<K,V>*)&other);
    }

    VariableSizeMap(void** data, unsigned i_capacity)
        : VariableSizeMapBase<K,V,growth,Hash,fill,_empty,_tombstone,false>(data, i_capacity)
    {
    }

    typedef VariableSizeMapEntry<K,V> Entry;

    inline void operator=(const VariableSizeMap<K,V>& other)
    {
        this->assign((VariableSizeMapBase<K,V>*)&other);
    }

    ~VariableSizeMap()
    {}
    
    
    inline bool tryGet(K key, V* o_value)
    {
        Entry* p = this->scan(key, false);
        if (p != NULL) {
            *o_value = p->value;
        }
        return p != NULL;
    }

    inline V* tryFind(K key)
    {
        Entry* p = this->scan(key, false);
        return p != NULL ? &p->value : NULL;
    }

    inline V get(K key)
    {
        Entry* p = this->scan(key, false);
        _ASSERT(p != NULL);
        return p->value;
    }

    bool erase(K key)
    {
        Entry* p = this->scan(key, false);
        if (p != NULL) {
            p->key = K(_tombstone);
            this->count--;
        }
        return p != NULL;
    }
    
    inline V& operator[](K key)
    {
        Entry* p = this->scan(key, false);
        _ASSERT(p != NULL);
        return p->value;
    }

    inline void put(K key, V value)
    {
        V* p;
        if (! tryAdd(key, value, &p)) {
            *p = value;
        }
    }

    inline V* getOrAdd(K key)
    {
        V* p = tryFind(key);
        if (p == NULL) {
            tryAdd(key, V(), &p);
        }
        return p;
    }

    inline bool tryAdd(K key, V value, V** o_pvalue)
    {
        while (true) {
            Entry* p = this->scan(key, true);
            if (p == NULL) {
                this->grow();
                p = this->scan(key, true);
                _ASSERT(p != NULL);
            }
            K prior = p->key;
            if (prior == key) {
                *o_pvalue = &p->value;
                return false;
            }
            if (prior == _empty || prior == _tombstone) {
                // single-threaded
                p->key = key;
                p->value = value;
                this->count++;
		        // hack!! to get around gcc bug
		        int o = this->occupied;
		        int l = this->limit;
                bool occupy = prior == _empty;
                if (o < l || ! occupy) {
                    this->occupied += occupy;
                    *o_pvalue = &p->value;
                } else {
                    this->grow();
                    *o_pvalue = &this->scan(key, false)->value; // lookup again after rehashing
                }
                return true;
            }
        }
    }

    void exchange(VariableSizeMap& other)
    {
        Entry* e = this->entries; this->entries = other.entries; other.entries = e;
        int x = this->capacity; this->capacity = other.capacity; other.capacity = x;
        x = this->count; this->count = other.count; other.count = x;
        x = this->limit; this->limit = other.limit; other.limit = x;
        x = this->occupied; this->occupied = other.occupied; other.occupied = x;
    }
};

template< typename K, typename V>
class VariableSizeMapBig
    : public VariableSizeMap<K,V,150,MapNumericHash<K>,80,0,-1,true>
{
public:
    VariableSizeMapBig(int n = 10000) : VariableSizeMap<K,V,150,MapNumericHash<K>,80,0,-1,true>(n) {}

    void assign(VariableSizeMapBig<K,V>* other)
    {
        // todo: avoid copying from base class, c++ inheritance is nonsensical
        if (this->entries != NULL) {
            BigDealloc(this->entries);
        }
        this->entries = other->entries;
        this->capacity = other->capacity;
        this->count = other->count;
        this->limit = other->limit;
        this->hash = other->hash;
        this->occupied = other->occupied;
        other->entries = NULL;
        other->count = 0;
    }
};

typedef VariableSizeMap<unsigned,unsigned> IdMap;
typedef VariableSizeMap<unsigned,int> IdIntMap;
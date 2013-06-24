#pragma once

#include "Compat.h"
#include "BigAlloc.h"
#include "VariableSizeVector.h"

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
    int fill = 90,
    int _empty = 0,
    int _tombstone = -1,
    int _busy = -2,
    bool multi = false>
class VariableSizeMapBase
{
protected:
    VariableSizeMapBase(int i_capacity = 16)
        : entries(NULL), count(0), capacity(i_capacity)
    {
        reserve(max(16,i_capacity));
    }

    VariableSizeMapBase(void** data, unsigned i_capacity)
        : entries((Entry*) (2 + (_int64*) *data)),
        capacity(i_capacity),
        count((int) ((_int64*)*data)[0]),
        limit((int) ((_int64*)*data)[1])
    {
        *data = ((char*)*data) + (size_t) i_capacity * sizeof(Entry) + 2 * sizeof(_int64);
    }
    
    inline void assign(VariableSizeMapBase<K,V>* other)
    {
        if (entries != NULL) {
            delete [] entries;
        }
        entries = other->entries;
        capacity = other->capacity;
        count = other->count;
        limit = other->limit;
        hash = other->hash;
        other->entries = NULL;
        other->count = 0;
    }

    inline void grow()
    {
        _int64 larger = ((_int64) capacity * growth) / 100;
        _ASSERT(larger < INT32_MAX);
        reserve((int) larger);
    }

public:
    inline int size()
    { return count; }

    inline int getCapacity()
    { return capacity; }
    
    ~VariableSizeMapBase()
    {
        if (entries != NULL) {
            delete [] entries;
        }
        entries = NULL;
        count = 0;
    }
    
    void reserve(int larger)
    {
        Entry* old = entries;
        int small = capacity;
        capacity = larger;
        entries = new Entry[larger];
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
                    Entry* p = scan(k, true);
                    _ASSERT(p != NULL);
                    p->key = k;
                    p->value = old[i].value;
                    count++;
                }
            }
            delete[] old;
        }
    }
    
    void clear()
    {
        if (entries != NULL) {
            if (_empty == 0) {
                // optimize zero case
                memset(entries, 0, capacity * sizeof(Entry));
            } else {
                const K e(_empty);
                for (int i = 0; i < capacity; i++) {
                    entries[i].key = e;
                }
            }
        }
        count = 0;
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
        Entry* p = scan(key, false);
        return p != NULL ? p : end();
    }

    void writeFile(LargeFileHandle* file)
    {
        _int64 x = (_int64) count;
        WriteLargeFile(file, &x, sizeof(_int64));
        x = (_int64) limit;
        WriteLargeFile(file, &x, sizeof(_int64));
        WriteLargeFile(file, entries, sizeof(Entry) * (size_t) capacity);
    }
    
protected:
    
    static const int MaxQuadraticProbes = 3;

    void init(int& pos, int& i, K key)
    {
        _ASSERT(key != _empty && key != _tombstone && key != _busy);
        pos = hash(key) % capacity;
        i = 1;
        if (entries == NULL) {
            entries = new Entry[capacity];
            clear();
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
            while (growth == 0 && p->key == _busy) {
                // spin
            }
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
    int limit; // current limit (capacity * fill / 100)
    Hash hash;
};

// 
// Single-valued map
//
template< typename K, typename V, int growth = 150, typename Hash = MapNumericHash<K>,
    int fill = 90, int _empty = 0, int _tombstone = -1, int _busy = -2 >
class VariableSizeMap
    : public VariableSizeMapBase<K,V,growth,Hash,fill,_empty,_tombstone,_busy,false>
{
public:
    VariableSizeMap(int i_capacity = 16)
        : VariableSizeMapBase<K,V,growth,Hash,fill,_empty,_tombstone,_busy,false>(i_capacity)
    {}

    VariableSizeMap(const VariableSizeMap<K,V>& other)
    {
        assign((VariableSizeMapBase<K,V>*)&other);
    }

    VariableSizeMap(void** data, unsigned i_capacity)
        : VariableSizeMapBase<K,V,growth,Hash,fill,_empty,_tombstone,_busy,false>(data, i_capacity)
    {
    }

    typedef VariableSizeMapEntry<K,V> Entry;

    inline void operator=(const VariableSizeMap<K,V>& other)
    {
        assign((VariableSizeMapBase<K,V>*)&other);
    }

    ~VariableSizeMap()
    {}
    
    
    inline bool tryGet(K key, V* o_value)
    {
        Entry* p = scan(key, false);
        if (p != NULL) {
            *o_value = p->value;
        }
        return p != NULL;
    }

    inline V* tryFind(K key)
    {
        Entry* p = scan(key, false);
        return p != NULL ? &p->value : NULL;
    }

    inline V get(K key)
    {
        Entry* p = scan(key, false);
        _ASSERT(p != NULL);
        return p->value;
    }

    bool erase(K key)
    {
        Entry* p = scan(key, false);
        if (p != NULL) {
            p->key = K(_tombstone);
            this->count--;
        }
        return p != NULL;
    }
    
    inline V& operator[](K key)
    {
        Entry* p = scan(key, false);
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
            Entry* p = scan(key, true);
            if (p == NULL) {
                this->grow();
                p = scan(key, true);
                _ASSERT(p != NULL);
            }
            K prior = p->key;
            if (prior == key) {
                *o_pvalue = &p->value;
                return false;
            }
            if (prior == _empty || prior == _tombstone) {
                if (growth != 0) {
                    // single-threaded
                    p->key = key;
                    p->value = value;
                    this->count++;
		    // hack!! to get around gcc bug
		    int c = this->count;
		    int l = this->limit;
                    if (c < l) {
                        *o_pvalue = &p->value;
                    } else {
                        this->grow();
                        *o_pvalue = &scan(key, false)->value; // lookup again after rehashing
                    }
                    return true;
                } else {
                    // multi-threaded
                    if (sizeof(K) == 4) {
#if 0
                        if (InterlockedCompareExchange32AndReturnOldValue((volatile _uint32*) &p->key, busy, prior) != prior) {
                            continue;
                        }
#endif
                    } else if (sizeof(K) == 8) {
#if 0
                        if (InterlockedCompareExchange64AndReturnOldValue((volatile _uint64*) &p->key, busy, prior) != prior) {
                            continue;
                        }
#endif
                    } else {
                        p->key = K(_busy);
                    }
                    p->value = value;
#if 0
                    InsertWriteBarrier();
#endif
                    p->key = key;
                    int incremented = (int) InterlockedIncrementAndReturnNewValue((volatile int*) &this->count);
                    _ASSERT(incremented <= limit);
                    *o_pvalue = &p->value;
                    return true;
                }
            }
        }
    }
};

typedef VariableSizeMap<unsigned,unsigned> IdMap;
typedef VariableSizeMap<unsigned,int> IdIntMap;
// 
// Single-valued map
//
template< typename K, typename V, int growth = 150, typename Hash = MapNumericHash<K>, int fill = 90, K empty = K(), K tombstone = K(-1), K busy = K(-2) >
class VariableSizeMultiMap
    : public VariableSizeMapBase<K,V,growth,Hash,fill,empty,tombstone,busy,true>
{
public:
    VariableSizeMultiMap(int i_capacity = 16)
        : VariableSizeMapBase<K,V,growth,Hash,fill,empty,tombstone,busy,true>(i_capacity)
    {}

    VariableSizeMultiMap(VariableSizeMultiMap& other)
        : VariableSizeMapBase<K,V,growth,Hash,fill,empty,tombstone,busy,true>(other.capacity)
    {
        this->assign(&other);
    }

    VariableSizeMultiMap(void** data, unsigned i_capacity)
        : VariableSizeMapBase<K,V,growth,Hash,fill,empty,tombstone,busy,true>(data, i_capacity)
    {
    }

    typedef VariableSizeMapEntry<K,V> Entry;

    inline void operator=(VariableSizeMultiMap<K,V> other)
    {
        this->assign(&other);
    }

    ~VariableSizeMultiMap()
    {}

    class valueIterator
    {
    public:
        bool hasValue()
        { return pos < map->capacity && map->entries[pos].key != empty; }

        Entry* operator*() const
        { _ASSERT(pos < map->capacity); return &map->entries[pos]; }
        
        Entry* operator->() const
        { _ASSERT(pos < map->capacity); return &map->entries[pos]; }
        
        void next()
        {
            if (hasValue()) {
                K k;
                do {
                    if (! map->advance(pos, i, key)) {
                        _ASSERT(false);
                    }
                } while ((k = map->entries[pos].key) != key && k != empty);
            }
        }

        valueIterator()
            : map(NULL), pos(0), i(0), key()
        {
        }

        valueIterator(const valueIterator& other)
            : map(other.map), pos(other.pos), i(other.i), key(other.key)
        {}

        void operator= (const valueIterator& other)
        {
            map = other.map;
            pos = other.pos;
            i = other.i;
            key = other.key;
        }

    private:
        valueIterator(VariableSizeMultiMap* i_map, K i_key)
            : map(i_map), key(i_key)
        {
            map->init(pos, i, key);
            K k = map->entries[pos].key;
            if (k != key && k != empty) {
                next(); // skip tombstones & other keys
            }
        }

        friend class VariableSizeMultiMap;

        VariableSizeMultiMap* map;
        int pos;
        int i;
        K key;
    };

    friend class valueIterator;
    
    inline valueIterator getAll(K key)
    {
        return valueIterator(this, key);
    }

    inline bool hasKey(K key)
    {
        return getAll(key).hasValue(); // todo: optimize
    }

    inline bool contains(K key, V value)
    {
        for (valueIterator i = getAll(key); i.hasValue(); i.next()) {
            if (i->value == value) {
                return true;
            }
        }
        return false;
    }

    // always add even if value exists for key
    inline void add(K key, V value)
    {
        if (this->count >= this->limit) {
            this->grow();
        }
        Entry* p = this->scan(key, true);
        if (p == NULL) {
            // full of tombstones
            this->grow();
            p = this->scan(key, true);
            _ASSERT(p != NULL);
        }
        p->key = key;
        p->value = value;
        this->count++;
    }

    // if key-value exists, return false; else add & return true
    inline bool put(K key, V value)
    {
        int pos; int i;
        init(pos, i, key);
        if (pos == this->capacity) {
            add(key, value);
            return true;
        }
        Entry* slot = NULL;
        while (pos != this->capacity) {
            Entry* p = &this->entries[pos];
            if (p->key == key) {
                if (p->value == value) {
                    return false;
                }
            } else if (p->key == tombstone && slot == NULL) {
                slot = p; // remember in case we don't find a match
            } else if (p->key == empty) {
                slot = p; // got to end with no match, use this slot
                break;
            }
            if (! this->advance(pos, i, key)) {
                break;
            }
        }
	    // hack!! to get around gcc bug
	    int c = this->count;
	    int l = this->limit;
        if (slot != NULL && c < l) {
            _ASSERT(slot->key == empty || slot->key == tombstone);
            slot->key = key;
            slot->value = value;
            this->count++;
        } else {
            this->add(key, value);
        }
        _ASSERT(this->contains(key, value));
        return true;
    }

    inline bool erase(K key, V value)
    {
        for (valueIterator i = getAll(key); i.hasValue(); i.next()) {
            if (i->value == value) {
                i->key = tombstone;
                this->count--;
                return true;
            }
        }
        return false;
    }
    
    inline int eraseAll(K key)
    {
        int n = 0;
        for (valueIterator i = getAll(key); i.hasValue(); i.next()) {
            i->key = tombstone;
            this->count--;
            n++;
        }
        return n;
    }

    // whether a's values are a subset of b's values
    inline bool isSubset(K a, K b)
    {
        for (valueIterator i = getAll(a); i.hasValue(); i.next()) {
            if (! contains(b, i->value)) {
                return false;
            }
        }
        return true;
    }

    inline bool intersects(K a, K b)
    {
        for (valueIterator i = getAll(a); i.hasValue(); i.next()) {
            if (contains(b, i->value)) {
                return true;
            }
        }
        return false;
    }

    inline void getAll(K key, VariableSizeVector<K>& o_result)
    {
        for (valueIterator i = getAll(key); i.hasValue(); i.next()) {
            o_result.push_back(i->value);
        }
    }
};

typedef VariableSizeMultiMap<unsigned,unsigned> IdMultiMap;

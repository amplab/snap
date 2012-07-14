#include "stdafx.h"
#include "Compat.h"

namespace sim {

//
// 'DiploidMutator' helper.
// 
class OffsetMap {
public:
    OffsetMap() : totalGain(0) {}
    //
    // Register a mutation at 'loc' that replaces 'before' bases with 'after' bases.
    // Mutations must be added in sorted order, and must not overlap.
    //
    void addMutation(unsigned loc, unsigned before, unsigned after) {
        _ASSERT(before);
        if (before != after) {
            unsigned replacementStart = loc + totalGain;
            unsigned replacementEnd = replacementStart + after;
            if (before > after)
                addEntry(replacementEnd, -(int)(before - after));
            else
                for (unsigned i = replacementStart + before; i < replacementEnd; ++i)
                    addEntry(i, 1);
        }
    }
    //
    // Coordinate conversion: mutated to unmutated.
    // 
    unsigned getTranslation(unsigned descendant) const {
        typedef Entries::const_reverse_iterator Iter;
        Iter previous = std::lower_bound(entries.rbegin(), entries.rend(),
                                         Entry(descendant, 42), greater);
        return descendant - (previous == entries.rend() ? 0 : previous->second);
    }
private:
    typedef std::pair<unsigned, int> Entry;
    typedef std::vector<Entry> Entries;
    static bool greater(const Entry &entry1, const Entry &entry2) {
        return entry1.first > entry2.first;
    }
    void addEntry(unsigned descendant, int localGain) {
        Entry newEntry(descendant, totalGain += localGain);
        _ASSERT(entries.empty() || greater(newEntry, entries.back()));
        entries.push_back(newEntry);
    }
    Entries entries;
    int totalGain;
};

}

/*++

Module Name:

    Cigar.h

Abstract:

    Code to parse cigar strings

Authors:

    Bill Bolosky, June, 2012

Environment:

    User mode service.

Revision History:


--*/

class CigarString {
public:
    //
    // Checks to see if the string is "*"
    //
    static bool IsCigarStringOmitted(const char *string);

    //
    // Can't call this on an omitted string.
    //
    CigarString(const char *string, size_t stringLength);

    unsigned        getLeftClip();
    unsigned        getRightClip();

    //
    // The tests for the state of the bases all use whichChar of the clipped read (ignoring the skipped portions).
    //

    bool            isBaseEqual(unsigned whichBase);
    bool            isBaseMismatch(unsigned whichBase); // This means SNP, not insert or delete
    unsigned        getCountOfDeletionsImmediatelyBeforeBase(unsigned whichBase);
    unsigned        getCountOfDeletionsAnywhereBeforeBase(unsigned whichBase);
    
private:

    //
    // The status of a base in a Read.  We don't have delete here, because if it's deleted then there's
    // no base.  Deletions are represented as coming between bases.  Similarly, there's no clipped.
    //
    enum BaseStatus {Equal, Mismatch, Insertion};

    unsigned readLength;
    BaseStatus *baseStatus;
    unsigned *totalDeletionsBefore; // This is the number of deletions anywhere in the read before this base
};
/*++

Module Name:

    Cigar.cpp

Abstract:

    Code to parse cigar strings

Authors:

    Bill Bolosky, June, 2012

Environment:

    User mode service.

Revision History:


--*/

#include "stdafx.h"
#include "Cigar.h"

    bool
CigarString::IsCigarStringOmitted(const char *string)
{
    return NULL == string || *string == '\0' || *string == '*';
}

CigarString::CigarString(const char *string, size_t stringLength)
{
    //
    // Start by counting the number of elements in the string.  An element is a count followed by a single character.
    //
    
}
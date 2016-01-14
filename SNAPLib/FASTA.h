/*++

Module Name:

    FASTA.h

Abstract:

    FASTA reader

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#pragma once

#include "Genome.h"

    const Genome *
ReadFASTAGenome(const char *fileName, const char *pieceNameTerminatorCharacters, bool spaceIsAPieceNameTerminator, unsigned chromosomePaddingSize, const char* chrTag, const char* chrMapFilename, AltContigMap* altMap);

//
// The FASTA appending functions return whether the write was successful.
// 
// WARNING: They write very long lines.
// According to Wikipedia, a FASTA file's line limit should be 120, or better, 79.
// Unix workaround if the piece names aren't too long: 'fold -w 79'.
//

    bool
AppendFASTAGenome(const Genome *, FILE *fasta);

// utility for parsing FASTA tags
    bool
FindFASTATagValue(const char* lineBuffer, const char* tag, const char ** pTagValue, int * pValueLength);

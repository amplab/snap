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

//
// There are several ways of specifying ALT contigs.  There is an opt-in list of ALTs, an opt-out list of regular chromosomes (these must be mutually
// exclusive), and a size cutoff below which is contig is an ALT.  The opt-in and opt-out lists supersede the size cutoff.
//

	const Genome *
ReadFASTAGenome(
	const char		*fileName,
	const char		*pieceNameTerminatorCharacters,
	bool			 spaceIsAPieceNameTerminator,
	unsigned		 chromosomePaddingSize,
	const char		*opt_in_alt_names,
	int				 opt_in_alt_names_count,
	const char		*opt_out_alt_names,
	int				 opt_out_alt_names_count,
	GenomeDistance	 maxSizeForAutomaticALT);

//
// The FASTA appending functions return whether the write was successful.
// 
// WARNING: They write very long lines.
// According to Wikipedia, a FASTA file's line limit should be 120, or better, 79.
// Unix workaround if the piece names aren't too long: 'fold -w 79'.
//

    bool
AppendFASTAGenome(const Genome *, FILE *fasta);

//
// This is arbitrary; is there some existing convention?
//
inline const char *diploidFASTASexPrefix(bool male)
{
    return male ? "PATERNAL|" : "MATERNAL|";
}

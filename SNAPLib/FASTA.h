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
	const char* const*opt_in_alt_names,
	int				 opt_in_alt_names_count,
	const char* const*opt_out_alt_names,
	int				 opt_out_alt_names_count,
	GenomeDistance	 maxSizeForAutomaticALT,
	bool             autoAlt,
    char            **alt_liftover_contig_names,
	unsigned		*alt_liftover_contig_flags,
	char			**alt_liftover_proj_contig_names,
	unsigned		*alt_liftover_proj_contig_offsets,
	char			**alt_liftover_proj_cigar,
	int				 alt_liftover_count);

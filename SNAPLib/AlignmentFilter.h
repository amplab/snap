/*++

Module Name:

    AlignmentFilter.h

Abstract:

    Filters transcriptome and genome alignments, allows for conversion from transcriptome
    to genome coordinates. Also handles fusion gene searching

Authors:

    Andrew Magis, June, 2013

Environment:

Revision History:

--*/

#pragma once

//System includes
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>

#include "stdafx.h"
#include "PairedEndAligner.h"
#include "BaseAligner.h"
#include "Genome.h"
#include "GTFReader.h"

using namespace std;

static const unsigned maxMAPQ = 70;

class Alignment {

    friend class AlignmentFilter;
    friend class AlignmentPair;

    public:
    
        Alignment(unsigned location_, Direction direction_, int score_, int mapq, string rname_, unsigned pos_, unsigned pos_end_, unsigned pos_original_, string transcript_id_, string gene_id_, bool isTranscriptome_);
        virtual ~Alignment() {}
        Alignment(const Alignment &rhs);
        Alignment& operator=(const Alignment&rhs);
        bool operator<(const Alignment &rhs) const;
        void Print();
        
    protected:
    
        unsigned location;
        Direction direction;
        int score;
        int mapq;
        string rname;
        unsigned pos;
        unsigned pos_end;
        unsigned pos_original;
        bool isTranscriptome;
        string transcript_id;
        string gene_id;
        string hashkey;
        
};

class AlignmentPair {

    friend class AlignmentFilter;

    public: 
    
        AlignmentPair(Alignment *align1_, Alignment *align2_, char flag_, bool is_unannotated_, bool is_backspliced_);
        virtual ~AlignmentPair() {}
        AlignmentPair(const AlignmentPair &rhs);
        AlignmentPair& operator=(const AlignmentPair &rhs);
        bool operator<(const AlignmentPair &rhs) const;
        void Print();
        
    protected:
    
        Alignment *align1;
        Alignment *align2;
        char flag;
        int distance;
        unsigned score;
        bool is_unannotated;
        bool is_backspliced;

};

typedef std::map<std::string, Alignment> alignment_map;

class AlignmentFilter {
  
    public:
    
        //Constructors/destructor
        AlignmentFilter(Read *read0_, Read *read1_, const Genome* genome_, const Genome* transcriptome_, GTFReader* gtf_, unsigned minSpacing_, unsigned maxSpacing_, unsigned confDiff_, unsigned maxDist_, unsigned seedLen_, BaseAligner *specialAligner, bool _enableFusions=false);
        virtual ~AlignmentFilter();  
        
        //Functions
        int AddAlignment(unsigned location, Direction direction, int score, int mapq, bool isTranscriptome, bool isMate0); 
        int Filter(PairedAlignmentResult* result);
        AlignmentResult FilterSingle(unsigned* location, Direction* direction, int* score, int* mapq, bool* isTranscriptome, unsigned* tlocation);   
        
    protected:
    
        void HashAlignment(Alignment& alignment, alignment_map& hashtable);
        void ConvertToGenomic();
        void ProcessPairs(PairedAlignmentResult* result, std::vector<AlignmentPair> &pairs);
        void CheckNoRC(PairedAlignmentResult* result, std::vector<AlignmentPair> &pairs);
        void FindPartialMatches(PairedAlignmentResult* result, AlignmentPair &pair);
        
        //Called on all unaligned reads to look for unknown splice junctions
        void UnalignedRead(alignment_map &mate, Read *read, unsigned minDiff);
        bool ProcessSplices(std::vector<AlignmentPair> &pairs, unsigned minDiff);
              
        //Temp printing
        void PrintMaps(seed_map& map, seed_map& mapRC);
        

        Read *read0;
        Read *read1;
        const Genome* genome;
        const Genome* transcriptome;
        GTFReader* gtf;
        alignment_map mate0;
        alignment_map mate1;        
        unsigned minSpacing;
        unsigned maxSpacing;
        unsigned confDiff;
        unsigned maxDist;
        unsigned seedLen;
        unsigned genome_mapq;
        BaseAligner *specialAligner;
        bool enableFusions;

};




#pragma once

//System includes
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>

#include "PairedEndAligner.h"
#include "Genome.h"
#include "GTFReader.h"

using namespace std;

class Alignment {

    friend class AlignmentFilter;

    public:
        Alignment(unsigned location_, bool isRC_, int score_, string rname_, unsigned pos_, unsigned pos_original_, string transcript_id_, bool isTranscriptome_);
        virtual ~Alignment();
        
        int Score() { return score; }
        void Print();
        
    protected:
    
        unsigned location;
        bool isRC;
        int score;
        string rname;
        int pos;
        int pos_original;
        bool isTranscriptome;
        string transcript_id;
        string hashkey;
        
};

typedef std::map<std::string, Alignment> alignment_map;
typedef std::pair<Alignment*, Alignment*> alignment_pair;

class AlignmentFilter {
  
    public:
    
        //Constructors/destructor
        AlignmentFilter(const Genome* genome_, const Genome* transcriptome_, GTFReader* gtf_, unsigned minSpacing_, unsigned maxSpacing_, unsigned confDiff_);
        virtual ~AlignmentFilter();  
        
        //Functions
        int AddAlignment(unsigned location, bool isRC, int score, bool isTranscriptome, bool isMate0); 
        int Filter(PairedAlignmentResult* result);    

    protected:
    
        int HashAlignment(Alignment& alignment, alignment_map& hashtable);
        void ConvertToGenomic();
                
        const Genome* genome;
        const Genome* transcriptome;
        GTFReader* gtf;
        alignment_map mate0;
        alignment_map mate1;
        unsigned minSpacing;
        unsigned maxSpacing;
        unsigned confDiff;


};

inline std::string ToString(const unsigned& arg)
{
  char buffer[65];
  sprintf(buffer, "%u", arg) ;
  return std::string(buffer);
}

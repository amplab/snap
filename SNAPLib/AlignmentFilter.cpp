#include "AlignmentFilter.h"

Alignment::Alignment(unsigned location_, bool isRC_, int score_, string rname_, unsigned pos_, unsigned pos_original_, string transcript_id_, bool isTranscriptome_) 
    : location(location_), isRC(isRC_), score(score_), rname(rname_), pos(pos_), pos_original(pos_original_), transcript_id(transcript_id_), isTranscriptome(isTranscriptome_)
{

    //Build the hashkey for this alignment
    hashkey = rname + '_' + ToString(pos);

}

Alignment::~Alignment() {}

void Alignment::Print() {
    printf("%u\t%d\t%d\t%s\t%u\t%s\t%d\n", location, isRC, score, rname.c_str(), pos, hashkey.c_str(), isTranscriptome);
}

AlignmentFilter::AlignmentFilter(const Genome* genome_, const Genome* transcriptome_, GTFReader* gtf_, unsigned minSpacing_, unsigned maxSpacing_, unsigned confDiff_) 
    : genome(genome_), transcriptome(transcriptome_), gtf(gtf_), minSpacing(minSpacing_), maxSpacing(maxSpacing_), confDiff(confDiff_)
{}

AlignmentFilter::~AlignmentFilter() {}

int AlignmentFilter::HashAlignment(Alignment& alignment, alignment_map& hashtable) {

    //Try to find this transcript in the transcript_map
    alignment_map::iterator pos = hashtable.find(alignment.hashkey);
        
    //If this sequence is not found, create a new vector to store this sequence (and others like it)
    if ((pos == hashtable.end())) {
        hashtable.insert(alignment_map::value_type(alignment.hashkey, alignment));
        return 0;
        
    //Place alignment with the better score
    } else {
        if (alignment.score < pos->second.score) {
            pos->second = alignment;
            return 1;
        }
    }
}

int AlignmentFilter::AddAlignment(unsigned location, bool isRC, int score, bool isTranscriptome, bool isMate0) {

    //Get the position and rname for this alignment
    string rname = "*";
    unsigned pos = 0;
    unsigned pos_original = 0;
    string transcript_id = "";
    
    //If this is, in fact, aligned to something
    if (location != 0xFFFFFFFF) {
    
        if (!isTranscriptome) {
    
            const Genome::Piece *piece = genome->getPieceAtLocation(location);
            rname = piece->name;
            pos_original = location - piece->beginningOffset + 1;
            pos = pos_original;

        //If we have a transcriptome read, convert the coordinates to genomic coordinates
        } else {
        
            const Genome::Piece *piece = transcriptome->getPieceAtLocation(location);
            rname = piece->name;           
            pos_original = location - piece->beginningOffset + 1;
            pos = pos_original;
                 
            //Convert the transcript rname and pos into genomic coordinates
            GTFTranscript& transcript = gtf->GetTranscript(rname);
            transcript_id = transcript.TranscriptID();
            rname = transcript.Chr();
            pos = transcript.GenomicPosition(pos);
            
        }
    }
    
    //Create the Alignment
    Alignment alignment(location, isRC, score, rname, pos, pos_original, transcript_id, isTranscriptome);

    //Add the alignment to the hash_table
    if (isMate0) {
        HashAlignment(alignment, mate0);           
    } else {
        HashAlignment(alignment, mate1);          
    }
    
    return 0;
}

bool AlignmentPairSort(const alignment_pair& pair1, const alignment_pair& pair2) {
    if (pair1.first->Score() + pair1.second->Score() < pair2.first->Score() + pair2.second->Score()) {
        return true;
    }
    return false;
}

int AlignmentFilter::Filter(PairedAlignmentResult* result) {

    std::vector<alignment_pair> pairs;
    
    /*
    printf("First\n");
    for (alignment_map::iterator m0 = mate0.begin(); m0 != mate0.end(); ++m0) {
        m0->second.Print();
    }
    printf("Second\n");
    for (alignment_map::iterator m1 = mate1.begin(); m1 != mate1.end(); ++m1) {
        m1->second.Print();
    }
    printf("\n");
    */
    
    //Iterate through all mate0;
    for (alignment_map::iterator m0 = mate0.begin(); m0 != mate0.end(); ++m0) {
    
        //Ensure this read is aligned
        if (m0->second.pos == 0) {
            continue;
        }
    
        //Iterate through all mate1
        for (alignment_map::iterator m1 = mate1.begin(); m1 != mate1.end(); ++m1) {
        
            //Ensure this read is aligned
            if (m1->second.pos == 0) {
                continue;
            }
            
            //Ensure both alignments are on the same chromosome
            if (m0->second.rname.compare(m1->second.rname) != 0) {
                //printf("failed1\n");
                //Possibility for fusion genes here
                continue;
            }
            
            //Ensure one alignment is reverse complemented
            if (((m0->second.isRC) && (m1->second.isRC)) ||
                (!m0->second.isRC) && (!m1->second.isRC)) {
                //printf("failed2\n");
                continue;
            }
            
            //This could be transcriptome or genomic distance
            //unsigned distance = abs(m0->second.pos_original-m1->second.pos_original);
            
            //This is the only valid comparison because this is absolute genomic distance.  
            //The other comparison could be way off depending on if we are comparing transcriptomic alignments
            //to genomic alignments, etc.  Just keep in mind that paired-end distances can be very great with 
            //spliced reads due to introns.  
            unsigned distance = abs(m0->second.pos-m1->second.pos);
            if ((distance < minSpacing) || (distance > maxSpacing)) {
                //printf("failed3: %d %d %d\n", distance, minSpacing, maxSpacing);
                continue;
            }
            
            //Create a paired-end mate
            pairs.push_back(alignment_pair(&m1->second, &m0->second));
            //pairs.push_back(alignment_pair(&m0->second, &m1->second));
        }
    }
    
    /*
    //Print the pairs
    printf("VALID PAIRS\n");
    for (vector<alignment_pair>::iterator it = pairs.begin(); it != pairs.end(); ++it) {
        it->first->Print();
        it->second->Print();
        printf("\n");
    }
    */

    //If there are no unique reads
    if (pairs.size() == 0) {
        
        //Paired end alignment not found
        result->status[0] = NotFound;
        result->status[1] = NotFound;    
        
        result->location[0] = 0;
        result->isRC[0] = 0;
        result->score[0] = 0;
        result->isTranscriptome[0] = false;
        
        result->location[1] = 0;
        result->isRC[1] = 0;
        result->score[1] = 0;
        result->isTranscriptome[1] = false;
        
        return 0;
    
    } else if (pairs.size() == 1) {
    
        /*
        printf("FINAL\n");
        pairs[0].first->Print();
        pairs[0].second->Print();
        printf("\n");
        */
    
        //Unique high quality hit
        result->status[0] = CertainHit;
        result->status[1] = CertainHit;

        result->location[0] = pairs[0].first->location;
        result->isRC[0] = pairs[0].first->isRC;
        result->score[0] = pairs[0].first->score;
        result->isTranscriptome[0] = pairs[0].first->isTranscriptome;
        
        result->location[1] = pairs[0].second->location;
        result->isRC[1] = pairs[0].second->isRC;
        result->score[1] = pairs[0].second->score;
        result->isTranscriptome[1] = pairs[0].second->isTranscriptome;
        
        return 1;
    
    } else {
         
        //Sort the scores by score
        sort(pairs.begin(), pairs.end(), AlignmentPairSort);
        
        /*
        printf("SORTED VALID PAIRS\n");
        for (vector<alignment_pair>::iterator it = pairs.begin(); it != pairs.end(); ++it) {
            it->first->Print();
            it->second->Print();
            printf("\n");
        }
        */
        
        result->location[0] = pairs[0].first->location;
        result->isRC[0] = pairs[0].first->isRC;
        result->score[0] = pairs[0].first->score;
        result->isTranscriptome[0] = pairs[0].first->isTranscriptome;
        
        result->location[1] = pairs[0].second->location;
        result->isRC[1] = pairs[0].second->isRC;
        result->score[1] = pairs[0].second->score;
        result->isTranscriptome[1] = pairs[0].second->isTranscriptome;
    
        //Check to see if the best alignment exceeds the second best alignment
        //by at least confDiff
        unsigned diff = (pairs[1].first->Score() + pairs[1].second->Score()) - 
                        (pairs[0].first->Score() + pairs[0].second->Score());
               
        if (diff >= confDiff) {

            //Unique high quality hit
            result->status[0] = CertainHit;
            result->status[1] = CertainHit;
            return 1;    
        
        } else {
        
            //Multiple hits
            result->status[0] = MultipleHits;
            result->status[1] = MultipleHits;
            return 2;    
        
        }  
    }
}


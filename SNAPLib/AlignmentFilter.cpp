/*++

Module Name:

    AlignmentFilter.cpp

Abstract:

    Filters transcriptome and genome alignments, allows for conversion from transcriptome
    to genome coordinates. Also handles fusion gene searching

Authors:

    Andrew Magis, June, 2013

Environment:

Revision History:

--*/

#include "stdafx.h"

#include "AlignmentFilter.h"

using std::min;

Alignment::Alignment(GenomeLocation location_, Direction direction_, int score_, int mapq_, string rname_, unsigned pos_, unsigned pos_end_, unsigned pos_original_, string transcript_id_, string gene_id_, bool isTranscriptome_) 
    : location(location_), direction(direction_), score(score_), mapq(mapq_), rname(rname_), pos(pos_), pos_end(pos_end_), pos_original(pos_original_), transcript_id(transcript_id_), gene_id(gene_id_), isTranscriptome(isTranscriptome_)
{

    //Build the hashkey for this alignment
    hashkey = rname + '_' + ToString(pos);

}

Alignment::Alignment(const Alignment &rhs)
    : location(rhs.location), direction(rhs.direction), score(rhs.score), mapq(rhs.mapq), rname(rhs.rname), pos(rhs.pos), pos_end(rhs.pos_end), pos_original(rhs.pos_original), transcript_id(rhs.transcript_id), gene_id(rhs.gene_id), isTranscriptome(rhs.isTranscriptome), hashkey(rhs.hashkey)    
{}

Alignment& Alignment::operator=(const Alignment &rhs) {
    if (this != &rhs) {
        location = rhs.location;
        direction = rhs.direction;
        score = rhs.score;
        mapq = rhs.mapq;
        rname = rhs.rname;
        pos = rhs.pos;
        pos_end = rhs.pos_end;
        pos_original = rhs.pos_original;
        isTranscriptome = rhs.isTranscriptome;
        transcript_id = rhs.transcript_id;
        gene_id = rhs.gene_id;
        hashkey = rhs.hashkey;
    }
    return *this;
}

bool Alignment::operator<(const Alignment &rhs) const {
    return (score < rhs.score);
}

void Alignment::Print() {
    printf("%u\t%d\t%d\t%d\t%s\t%u\t%s\t%d\n", location, direction, score, mapq, rname.c_str(), pos, hashkey.c_str(), isTranscriptome);
}

AlignmentPair::AlignmentPair(Alignment *align1_, Alignment *align2_, char flag_, bool is_unannotated_, bool is_backspliced_) 
    : align1(align1_), align2(align2_), flag(flag_), distance(0), is_unannotated(is_unannotated_), is_backspliced(is_backspliced_)
{

    //Calculate the score
    score = align1->score + align2->score;

    //If alignment1 is reverse complemented and alignment2 is not
    if ((align1->direction) && (!align2->direction)) {
        distance = align1->pos - align2->pos;
    } else if ((!align1->direction) && (align2->direction)) {
        distance = align2->pos - align1->pos;
    }
    
}    
    
AlignmentPair::AlignmentPair(const AlignmentPair &rhs) 
    : align1(rhs.align1), align2(rhs.align2), flag(rhs.flag), distance(rhs.distance), score(rhs.score), is_unannotated(rhs.is_unannotated), is_backspliced(rhs.is_backspliced)
{}

AlignmentPair& AlignmentPair::operator=(const AlignmentPair &rhs) {
    if (this != &rhs) {
        align1 = rhs.align1;
        align2 = rhs.align2;
        flag = rhs.flag;
        distance = rhs.distance;
        score = rhs.score;
        is_unannotated = rhs.is_unannotated;
        is_backspliced = rhs.is_backspliced;
    }
    return *this;
}

bool AlignmentPair::operator<(const AlignmentPair &rhs) const {
    return (score < rhs.score);
}

void AlignmentPair::Print() {
    printf("Alignment distance: %d %d %d\n", distance, align1->direction, align2->direction);
    align1->Print();
    align2->Print();
    printf("\n");
}

AlignmentFilter::AlignmentFilter(Read *read0_, Read *read1_, const Genome* genome_, const Genome* transcriptome_, GTFReader* gtf_, unsigned minSpacing_, unsigned maxSpacing_, unsigned confDiff_, unsigned maxDist_, unsigned seedLen_, BaseAligner *specialAligner_, bool _enableFusions) 
    : read0(read0_), read1(read1_), genome(genome_), transcriptome(transcriptome_), gtf(gtf_), minSpacing(minSpacing_), maxSpacing(maxSpacing_), confDiff(confDiff_), maxDist(maxDist_), seedLen(seedLen_), specialAligner(specialAligner_), genome_mapq(maxMAPQ), enableFusions(_enableFusions)
{}

AlignmentFilter::~AlignmentFilter() {}

void AlignmentFilter::HashAlignment(Alignment& alignment, alignment_map& hashtable) {

    //Try to find this transcript in the transcript_map
    alignment_map::iterator pos = hashtable.find(alignment.hashkey);
        
    //If this sequence is not found, create a new vector to store this sequence (and others like it)
    if ((pos == hashtable.end())) {
        hashtable.insert(alignment_map::value_type(alignment.hashkey, alignment));
        
    } else {
    
        //Place alignment with the better score
        if (alignment.score < pos->second.score) {
            pos->second = alignment;
        
        //If the two scores are equal, keep the transcriptome alignment (if any)
        } else if (alignment.score == pos->second.score) {
            if (alignment.isTranscriptome) {
                pos->second = alignment;
            }
        }
    }
}

int AlignmentFilter::AddAlignment(GenomeLocation location, Direction direction, int score, int mapq, bool isTranscriptome, bool isMate0) {

    //Get the position and rname for this alignment
    string rname = "*";
    unsigned pos = 0;
    unsigned pos_end = 0;
    unsigned pos_original = 0;
    string transcript_id = "";
    string gene_id = "";
    
    //If this is, in fact, aligned to something
    if (location != 0xFFFFFFFF) {
    
        if (!isTranscriptome) {
    
            const Genome::Contig *piece = genome->getContigAtLocation(location);

            //Check that piece is not NULL
            if (piece != NULL) {

                rname = piece->name;
                pos_original = location - piece->beginningLocation + 1;
                pos = pos_original;
            
                if (isMate0) {
                    pos_end = pos+read1->getDataLength()-1;
                } else {
                    pos_end = pos+read0->getDataLength()-1;
                }
            }
        
        //If we have a transcriptome read, convert the coordinates to genomic coordinates
        } else {
        
            const Genome::Contig *piece = transcriptome->getContigAtLocation(location);

            //Check that piece is not NULL
            if (piece != NULL) {
            
                rname = piece->name;           
                pos_original = location - piece->beginningLocation + 1;
                pos = pos_original;
            
                //Convert the transcript rname and pos into genomic coordinates
                const GTFTranscript& transcript = gtf->GetTranscript(rname);
                transcript_id = transcript.TranscriptID();
                gene_id = transcript.GeneID();
                        
                //RName is the chromosome name in genomic coordinates!
                rname = transcript.Chr();
            
                if (isMate0) {
                    pos_end = transcript.GenomicPosition(pos+read1->getDataLength()-1, 0);
                    pos = transcript.GenomicPosition(pos, read1->getDataLength());     
                } else {
                    pos_end = transcript.GenomicPosition(pos+read0->getDataLength()-1, 0);
                    pos = transcript.GenomicPosition(pos, read0->getDataLength()); 
                }
            }
        }
    } 
    
    //If the genomic location is valid
    if (pos != 0) {
    
        //Create the Alignment
        Alignment alignment(location, direction, score, mapq, rname, pos, pos_end, pos_original, transcript_id, gene_id, isTranscriptome);
        
        //Add the alignment to the hash_table
        if (isMate0) {
            HashAlignment(alignment, mate0);           
        } else {
            HashAlignment(alignment, mate1);          
        }
        
    }
    
    return 0;
}

AlignmentResult AlignmentFilter::FilterSingle(GenomeLocation* location, Direction* direction, int* score, int* mapq, bool* isTranscriptome, unsigned* tlocation) {

    std::vector<Alignment> alignments;

    //Here we simply pick the best alignment and return it
    for (alignment_map::iterator m0 = mate0.begin(); m0 != mate0.end(); ++m0) {
    
        //If this read passes the maxDist cutoff
        if (m0->second.score > (int) maxDist) {
            continue;
        }
    
        alignments.push_back(m0->second);
    }
    
    //If there are no alignments, then return NotFound
    if (alignments.size() == 0) {
        *location = 0;
        *direction = FORWARD;
        *score = 0;
        *mapq = 0;
        *isTranscriptome = false;
        *tlocation = 0;
        return NotFound;
        
    } else if (alignments.size() == 1) {
 
		if (alignments[0].isTranscriptome) {
			*tlocation = alignments[0].location;
			GenomeLocation offset;
			genome->getLocationOfContig(alignments[0].rname.c_str(), &offset);
			offset += alignments[0].pos - 1;
			*location = offset;

		} else {
			*location = alignments[0].location;
			*tlocation = 0;
		}

        *direction = alignments[0].direction;
        *score = alignments[0].score;
        *mapq = min(maxMAPQ, genome_mapq);
        *isTranscriptome = alignments[0].isTranscriptome;

        if (alignments[0].isTranscriptome) {
          gtf->IncrementReadCount(alignments[0].transcript_id, alignments[0].pos_original, alignments[0].pos, read1->getDataLength());
        }
        return SingleHit;
        
    } else {
        
        //Sort the alignments
        std::sort(alignments.begin(), alignments.end());
        
		if (alignments[0].isTranscriptome) {
			*tlocation = alignments[0].location;
			GenomeLocation offset;
			genome->getLocationOfContig(alignments[0].rname.c_str(), &offset);
			offset += alignments[0].pos - 1;
			*location = offset;
		} else {
			*location = alignments[0].location;
			*tlocation = 0;
		}
        
        *direction = alignments[0].direction;
        *score = alignments[0].score;
        //*mapq = alignments[0].mapq;
        *isTranscriptome = alignments[0].isTranscriptome;

        unsigned diff = alignments[1].score - alignments[0].score;              
        if (diff >= confDiff) {
            *mapq = min(maxMAPQ, genome_mapq);
  
            if (alignments[0].isTranscriptome) {
               gtf->IncrementReadCount(alignments[0].transcript_id, alignments[0].pos_original, alignments[0].pos, read1->getDataLength());
            }
            return SingleHit; 
        
        } else {
            *mapq = 1;
            return MultipleHits;    
        }
    }
}

int AlignmentFilter::Filter(PairedAlignmentResult* result) {

    std::vector<AlignmentPair> no_rc;
    std::vector<AlignmentPair> intragene_pairs;
    std::vector<AlignmentPair> intrachromosomal_pairs;
    std::vector<AlignmentPair> interchromosomal_pairs;
        
    /*
    printf("Align1\n");
    for (alignment_map::iterator m0 = mate0.begin(); m0 != mate0.end(); ++m0) {
        m0->second.Print();
    }
    printf("Align2\n");
    for (alignment_map::iterator m1 = mate1.begin(); m1 != mate1.end(); ++m1) {
        m1->second.Print();
    } 
    */
   
    char flag = 0;
    unsigned best_score = 10000;
    
    if ((mate0.size() == 0) && (mate1.size() == 0)) {
        flag |= 1 << FIRST_NOT_ALIGNED;
        flag |= 1 << SECOND_NOT_ALIGNED;
        
        //You ought to eliminate this and instead only allow single matches to dictate splicing

        if (enableFusions) {
            //UnalignedRead(read1, seedLen);
            //UnalignedRead(read0, seedLen);
        }

    //If there are no alignments for mate0
    } else if (mate0.size() == 0) {
        
        flag |= 1 << FIRST_NOT_ALIGNED;

        if (enableFusions) {
            UnalignedRead(mate1, read1, seedLen);
        }    

    //If there are no alignments for mate1
    } else if (mate1.size() == 0) {
    
        flag |= 1 << SECOND_NOT_ALIGNED;
     
        if (enableFusions) {
            UnalignedRead(mate0, read0, seedLen);    
        }
  
    }
    
    //Iterate through all mate0;
    for (alignment_map::iterator m0 = mate0.begin(); m0 != mate0.end(); ++m0) {
     
        if (m0->second.score > (int) maxDist) {
            continue;
        }
       
        //Iterate through all mate1
        for (alignment_map::iterator m1 = mate1.begin(); m1 != mate1.end(); ++m1) {
        
            //If this read passes the maxDist cutoff
            if (m1->second.score > (int) maxDist) {
                continue;
            }
            
            //Calculate distances between these two alignments
            int distance = 0;
            if ((m0->second.direction) && (!m1->second.direction)) {
                distance = m0->second.pos - m1->second.pos;
            } else if ((!m0->second.direction) && (m1->second.direction)) {
                distance = m1->second.pos - m0->second.pos;
            }
            
            //Check to see if this pair is 'circularized'
            bool is_backspliced = false;
            if (distance < -100) {
                is_backspliced = true;
            }
                  
            //Set the initial flag  
            flag = 0;
             
            //Ensure one alignment is reverse complemented
            if (((m0->second.direction) && (m1->second.direction)) ||
                ((!m0->second.direction) && (!m1->second.direction))) {
                flag |= 1 << NOT_REVERSE_COMPLIMENTED;
                no_rc.push_back(AlignmentPair(&m1->second, &m0->second, flag, false, is_backspliced));
                continue;
                           
            //If both reads are aligned to the transcriptome
            } else if ((m0->second.isTranscriptome) && (m1->second.isTranscriptome)) {
            
                //If they are on different chromosomes
                if (m0->second.rname.compare(m1->second.rname) != 0) {
                    flag |= 1 << ALIGNED_DIFF_CHR;
                    interchromosomal_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, false, is_backspliced));
                    continue;              
                    
                //If they are on the same chromosome and within within the gene boundary
                } else if (gtf->GetGene(m0->second.gene_id).CheckBoundary(m1->second.rname, m1->second.pos)) {
                    flag |= 1 << ALIGNED_SAME_GENE; 
                    intragene_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, false, is_backspliced));  
                    continue;
                    
                //If they are on the same chromosome, not within the same gene, but within the gene boundary
                } else if (gtf->GetGene(m1->second.gene_id).CheckBoundary(m0->second.rname, m0->second.pos)) {
                    flag |= 1 << ALIGNED_SAME_GENE; 
                    intragene_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, false, is_backspliced));  
                    continue;
                    
                //If they are on the same chromosome, not within the same gene, and not within the gene boundary, 
                } else {
                    flag |= 1 << ALIGNED_SAME_CHR;
                    intragene_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, false, is_backspliced));
                    continue;
                }    
            
            //If only mate1 is aligned to transcriptome
            } else if (m0->second.isTranscriptome) {
            
                //If they are on different chromosomes
                if (m0->second.rname.compare(m1->second.rname) != 0) {
 
                    flag |= 1 << ALIGNED_DIFF_CHR;
                    interchromosomal_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, true, is_backspliced));
                    continue;               

                //If they are on the same chromosome and within the gene boundary
                } else if (gtf->GetGene(m0->second.gene_id).CheckBoundary(m1->second.rname, m1->second.pos)) {
                
                    flag |= 1 << ALIGNED_SAME_GENE; 
                    intragene_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, true, is_backspliced));  
                    continue;
                
                //If they are on the same chromosome but not within the gene boundary
                } else {
                
                    flag |= 1 << ALIGNED_SAME_CHR; 
                    intragene_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, true, is_backspliced));
                    continue;                
                
                }
                
            //If only mate2 is aligned to transcriptome
            } else if (m1->second.isTranscriptome) {

                //If they are on different chromosomes
                if (m0->second.rname.compare(m1->second.rname) != 0) {
 
                    flag |= 1 << ALIGNED_DIFF_CHR;
                    interchromosomal_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, true, is_backspliced));
                    continue;     

                //If they are on the same chromosome and within the gene boundary
                } else if (gtf->GetGene(m1->second.gene_id).CheckBoundary(m0->second.rname, m0->second.pos)) {
                
                    flag |= 1 << ALIGNED_SAME_GENE; 
                    intragene_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, true, is_backspliced));
                    continue;
                    
                //If they are on the same chromosome but not within the gene boundary
                } else {
  
                    flag |= 1 << ALIGNED_SAME_CHR; 
                    intragene_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, true, is_backspliced));
                    continue;                   
                
                }
            
            //If neither are aligned to transcriptome, we can't be sure
            } else {
		
		
                //TESTME: TEST EFFECT OF THIS LINE INSTEAD OF THE CODE BLOCK BELOW
                //intragene_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, true, is_backspliced));
 
                
                //If they are on different chromosomes
                if (m0->second.rname.compare(m1->second.rname) != 0) {
 
                    flag |= 1 << ALIGNED_DIFF_CHR;
                    interchromosomal_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, true, is_backspliced));
                    continue;    
            
                } else {
                        
                    //Query the GTF interval tree for all genes overlapping this position
                    std::vector<GTFGene> results;
                    gtf->IntervalGenes(m0->second.rname, m0->second.pos, m0->second.pos, results);
  
                    //For each gene found, look within gene boundary for other read
                    bool found = false;
                    for (std::vector<GTFGene>::iterator it = results.begin(); it != results.end(); ++it) {
                        if (gtf->GetGene(it->GeneID()).CheckBoundary(m1->second.rname, m1->second.pos)) {
                            flag |= 1 << ALIGNED_SAME_GENE;
                            intragene_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, true, is_backspliced)); 
                            found = true; 
                            break;
                        }
                    }   
  
                    if (found) {
                        continue;
                    } else {
  
                        flag |= 1 << ALIGNED_SAME_CHR;
                        intragene_pairs.push_back(AlignmentPair(&m1->second, &m0->second, flag, true, is_backspliced));
                        continue;    
                    }        
                } 
            }
        }
    }
    
    /*
    printf("No RC Pairs: %u\n", no_rc.size());
    for (vector<AlignmentPair>::iterator it = no_rc.begin(); it != no_rc.end(); ++it) {
        it->Print();
    }
    
    printf("Same Gene Pairs: %u\n", intragene_pairs.size());
    for (vector<AlignmentPair>::iterator it = intragene_pairs.begin(); it != intragene_pairs.end(); ++it) {
        it->Print();
    }

    printf("Same Chr Pairs\n");
    for (vector<AlignmentPair>::iterator it = intrachromosomal_pairs.begin(); it != intrachromosomal_pairs.end(); ++it) {
        it->Print();
    }
    
    printf("Diff Chr Pairs\n");
    for (vector<AlignmentPair>::iterator it = interchromosomal_pairs.begin(); it != interchromosomal_pairs.end(); ++it) {
        it->Print();
    }
    */

    //Gene pairs always get priority.  If there is a paired end alignment
    if (intragene_pairs.size() > 0) {
    
        //Determine if these alignments are unique or not.
        ProcessPairs(result, intragene_pairs);
        
        /*
        //If this is still a good hit, check to make sure there is no partial hit that is better
        if (result->status[0] == SingleHit) {
            FindPartialMatches(result, no_rc[0]);
        }
        */

        //Here we check for negative reads which indicate circular RNAs
        if (result->status[0] == SingleHit) {
        
            /*
            //Here we check for circularized pairs, as well as pairs that map to unannotated regions            
             if (intragene_pairs[0].is_unannotated) {
                 result->flag[0] |= 1 << UNANNOTATED;
                 result->flag[1] |= 1 << UNANNOTATED;
                 gtf->IntrageneUnannotatedPair(intragene_pairs[0].align1->rname, 
                                               intragene_pairs[0].align1->pos, 
                                               intragene_pairs[0].align1->pos_end,
                                               intragene_pairs[0].align2->rname, 
                                               intragene_pairs[0].align2->pos, 
                                               intragene_pairs[0].align2->pos_end, 
                                               string(read0->getId(), read0->getIdLength()));              
             }
             
             if (intragene_pairs[0].is_backspliced) {
             
                 result->flag[0] |= 1 << CIRCULAR;
                 result->flag[1] |= 1 << CIRCULAR;
                 gtf->IntrageneCircularPair(intragene_pairs[0].align1->rname, 
                                            intragene_pairs[0].align1->pos, 
                                            intragene_pairs[0].align1->pos_end,
                                            intragene_pairs[0].align2->rname, 
                                            intragene_pairs[0].align2->pos, 
                                            intragene_pairs[0].align2->pos_end, 
                                            string(read0->getId(), read0->getIdLength()));                 
             
             }
             */
            
            //Add as count to GTF as well
            gtf->IncrementReadCount(intragene_pairs[0].align1->transcript_id, 
                                    intragene_pairs[0].align1->pos_original, 
                                    intragene_pairs[0].align1->pos, 
                                    read1->getDataLength(),
                                    intragene_pairs[0].align2->transcript_id, 
                                    intragene_pairs[0].align2->pos_original, 
                                    intragene_pairs[0].align2->pos, 
                                    read0->getDataLength());
        }
        
        result->fromAlignTogether = false;
        result->alignedAsPair = true;
        //result->nanosInAlignTogether;
        //result->nLVCalls;
        //result->nSmallHits;
        return 1;
    }
        
    //Pairs on the same chromosome get next priority
    if (intrachromosomal_pairs.size() > 0) {

        ProcessPairs(result, intrachromosomal_pairs);
    
        //If this is a good hit, check to make sure there is no RC hit that is better
        //if (result->status[0] == SingleHit) {
        //    CheckNoRC(result, no_rc);
        //}
        
        //If this pair is within some reasonable distance, then allow it
        //if (intrachromosomal_pairs[0].distance <= maxSpacing) {
        //    return 1;
        //}
        
        //If this is still a good hit, check to make sure there is no partial hit that is better
        if (result->status[0] == SingleHit) { 
            FindPartialMatches(result, intrachromosomal_pairs[0]);
        }

        //If this is still a good hit, add this in as a chr link
        if (result->status[0] == SingleHit) {

            //Link these positions in the GTF object
            gtf->IntrachromosomalPair(intrachromosomal_pairs[0].align1->rname, 
                                      intrachromosomal_pairs[0].align1->pos, 
                                      intrachromosomal_pairs[0].align1->pos_end,
                                      intrachromosomal_pairs[0].align2->rname, 
                                      intrachromosomal_pairs[0].align2->pos, 
                                      intrachromosomal_pairs[0].align2->pos_end, 
                                      string(read0->getId(), read0->getIdLength()), 
                                      string(read0->getData(), read0->getDataLength()),
                                      string(read1->getData(), read1->getDataLength()));      
        }
        
        result->fromAlignTogether = false;
        result->alignedAsPair = false;
        //result->nanosInAlignTogether;
        //result->nLVCalls;
        //result->nSmallHits;

        return 1;
    }    

    //Next comes different chromosome pairs
    if (interchromosomal_pairs.size() > 0) {
    
        ProcessPairs(result, interchromosomal_pairs);
        
        //If this is a good hit, check to make sure there is no RC hit that is better
        //if (result->status[0] == SingleHit) {
        //    CheckNoRC(result, no_rc);
        //}
       
        //If this is still a good hit, check to make sure there is no partial hit that is better
        if (result->status[0] == SingleHit) { 
            FindPartialMatches(result, interchromosomal_pairs[0]);
        }
        
        //If this is still a good hit, add this in as a gene link
        if (result->status[0] == SingleHit) {

            //Link these positions in the GTF object
            gtf->InterchromosomalPair(interchromosomal_pairs[0].align1->rname, 
                                      interchromosomal_pairs[0].align1->pos, 
                                      interchromosomal_pairs[0].align1->pos_end,
                                      interchromosomal_pairs[0].align2->rname, 
                                      interchromosomal_pairs[0].align2->pos, 
                                      interchromosomal_pairs[0].align2->pos_end, 
                                      string(read0->getId(), read0->getIdLength()), 
                                      string(read0->getData(), read0->getDataLength()),
                                      string(read1->getData(), read1->getDataLength()));
        }
        
        result->fromAlignTogether = false;
        result->alignedAsPair = false;
        //result->nanosInAlignTogether;
        //result->nLVCalls;
        //result->nSmallHits;
        
        return 1;
    }
   
    //Finally, we check the noRC alignments.  This is to pick up any fusion genes that 
    //have switched orientations relative to each other
    if (no_rc.size() > 0) {

        ProcessPairs(result, no_rc);
     
        //If this is still a good hit, check to make sure there is no partial hit that is better
        if (result->status[0] == SingleHit) {
            FindPartialMatches(result, no_rc[0]);
        }

        //If this is still a good hit, add this in as a gene link
        if (result->status[0] == SingleHit) {

            if (no_rc[0].align1->rname.compare(no_rc[0].align2->rname) == 0) {
        
              //Link these positions in the GTF object
              gtf->IntrachromosomalPair(no_rc[0].align1->rname, 
                                        no_rc[0].align1->pos, 
                                        no_rc[0].align1->pos_end,
                                        no_rc[0].align2->rname, 
                                        no_rc[0].align2->pos, 
                                        no_rc[0].align2->pos_end,
                                        string(read0->getId(), read0->getIdLength()), 
                                        string(read0->getData(), read0->getDataLength()),
                                        string(read1->getData(), read1->getDataLength()));

            } else {

                 //Link these positions in the GTF object              
                 gtf->InterchromosomalPair(no_rc[0].align1->rname, 
                                           no_rc[0].align1->pos, 
                                           no_rc[0].align1->pos_end, 
                                           no_rc[0].align2->rname, 
                                           no_rc[0].align2->pos, 
                                           no_rc[0].align2->pos_end,
                                           string(read0->getId(), read0->getIdLength()), 
                                           string(read0->getData(), read0->getDataLength()),
                                           string(read1->getData(), read1->getDataLength()));
            }
        }

        result->fromAlignTogether = false;
        result->alignedAsPair = false;
        //result->nanosInAlignTogether;
        //result->nLVCalls;
        //result->nSmallHits;

        return 1;
    }
        
    //If we do not have any gene pairs or nogene pairs, for now we will 
    //return a bad alignment

    result->tlocation[0] = 0;
    result->status[0] = NotFound;
    result->location[0] = 0;
    result->direction[0] = FORWARD;
    result->score[0] = 0;
    result->mapq[0] = 0;
    result->isTranscriptome[0] = false;
    
    result->tlocation[1] = 0;
    result->status[1] = NotFound;  
    result->location[1] = 0;
    result->direction[1] = FORWARD;
    result->score[1] = 0;
    result->mapq[1] = 0;
    result->isTranscriptome[1] = false;
    
    result->fromAlignTogether = false;
    result->alignedAsPair = false;
    //result->nanosInAlignTogether;
    //result->nLVCalls;
    //result->nSmallHits;
    
    return 0;    
          
}

void AlignmentFilter::UnalignedRead(alignment_map &mate, Read *read, unsigned minDiff) {

    seed_map map, mapRC;
    
    //Vectors to store potential splices
    std::vector<AlignmentPair> intragene_unannotated_splices;
    std::vector<AlignmentPair> intrachromosomal_splices;
    std::vector<AlignmentPair> interchromosomal_splices;
    
    specialAligner->setReadId(0);
    specialAligner->CharacterizeSeeds(read, map, mapRC);   
    char flag = 0;
    
    //PrintMaps(map, mapRC);

	std::vector<Alignment*> alignments;
        bool found = false;

        for (seed_map::iterator it = map.begin(); it != map.end(); ++it) {

                //Convert both segments to genomic coordinates
                const Genome::Contig *piece0 = genome->getContigAtLocation(it->first);
                string chr0 = piece0->name;
                int pos0 = it->first - piece0->beginningLocation + 1;

                //Make sure this location is near to at least one mate
                for (alignment_map::iterator mit = mate.begin(); mit != mate.end(); ++mit) {

                    //printf("%s %s %u %u\n", mit->second.rname.c_str(), chr0.c_str(), mit->second.pos, pos0);

                    if (((mit->second.rname.compare(chr0)) == 0) &&
                        (abs((int)mit->second.pos - (int)pos0) < 1000000)) {
                        found = true;
                        break;
                    }
                }
          }

        if (!found) {
        for (seed_map::iterator it = mapRC.begin(); it != mapRC.end(); ++it) {

                //Convert both segments to genomic coordinates
                const Genome::Contig *piece0 = genome->getContigAtLocation(it->first);
                string chr0 = piece0->name;
                int pos0 = it->first - piece0->beginningLocation + 1;

                //Make sure this location is near to at least one mate
                for (alignment_map::iterator mit = mate.begin(); mit != mate.end(); ++mit) {

                    //printf("%s %s %u %u\n", mit->second.rname.c_str(), chr0.c_str(), mit->second.pos, pos0);

                    if (((mit->second.rname.compare(chr0)) == 0) &&
                        (abs((int)mit->second.pos - (int)pos0) < 1000000)) {
                        found = true;
                        break;
                    }
                }
            }
       }

       if (found) {

	//Create partial alignments from map
	for (seed_map::iterator it = map.begin(); it != map.end(); ++it) {

		// Get length of mapping segment
		unsigned length0 = (*(it->second).rbegin() - *(it->second).begin()) + seedLen;

		//Convert segment to genomic coordinates
		const Genome::Contig *piece0 = genome->getContigAtLocation(it->first);
		string chr0 = piece0->name;
		int pos0 = it->first - piece0->beginningLocation + 1; 
	
		//Calculate the consecutive region of the genome that contains this segment
		unsigned start0 = pos0 + *(it->second).begin();
		unsigned end0 = start0 + length0 - 1;
		//printf("0 [%s:%u-%u]\n", chr0.c_str(), start0, end0);
										   
		//Create new alignments for each segment
                //printf("%s:%u-%u [%u]\n", chr0.c_str(), start0, end0, length0);
		alignments.push_back(new Alignment(it->first, false, length0, 0, chr0, start0, end0, start0, "transcript_id", "gene_id", false));
		
	}

	//Create partial alignments from mapRC
	for (seed_map::iterator it = mapRC.begin(); it != mapRC.end(); ++it) {
			
		//Now we check to see if the two lengths are nearly the length of the read
		unsigned length0 = (*(it->second).rbegin() - *(it->second).begin()) + seedLen;
						
		//Convert both segments to genomic coordinates
		const Genome::Contig *piece0 = genome->getContigAtLocation(it->first);
		string chr0 = piece0->name;
		int pos0 = it->first - piece0->beginningLocation + 1; 

		//Calculate the consecutive region of the genome that contains this segment
		unsigned start0 = pos0 + read->getDataLength() - (*(it->second).rbegin() + seedLen);
		unsigned end0 = start0 + length0 - 1;
		//printf("RC0 [%s:%u-%u]\n", chr0.c_str(), start0, end0);
		
		//Create new alignments for each segment
                //printf("%s:%u-%u [%u]\n", chr0.c_str(), start0, end0, length0);
		alignments.push_back(new Alignment(it->first, true, length0, 0, chr0, start0, end0, start0, "transcript_id", "gene_id", false));

	}
        }

	//Now we go through all partial alignments that have been created
	for (std::vector<Alignment*>::iterator it0 = alignments.begin(); it0 != alignments.end(); ++it0) {
		for (std::vector<Alignment*>::iterator it1 = it0; it1 != alignments.end(); ++it1) {
			
	    Alignment *align0 = (*it0);
            Alignment *align1 = (*it1);
	
            //Do not compare same alignments
            if (align0 == align1) {
                continue;
            }	
            
            //printf("SCORES: %u %u %u %u\n", align0->score, align1->score, read->getDataLength(), seedLen);

            //If not enough of the read is represented
            if ((align0->score + align1->score) < (int) (read->getDataLength() - seedLen)) {
                continue;
            }		
            
            bool is_backspliced = false;
            //If one segment begins after the other one ends
            if (align0->pos > align1->pos_end) {
            	is_backspliced = true;  //this is not actually correct anymore
            } else if (align1->pos > align0->pos_end) {
            	is_backspliced = false;  //this is not actually correct anymore
            } else {
            	continue;
            }	  
            
            //If they are on different chromosomes
            if (align0->rname.compare(align1->rname) != 0) {
            
                interchromosomal_splices.push_back(AlignmentPair(align0, align1, flag, true, is_backspliced));
                continue;    
        
            } else {
            
                //Query the GTF interval tree for all genes overlapping this position
                std::vector<GTFGene> results;
                gtf->IntervalGenes(align0->rname, align0->pos, align0->pos_end, results);

                //For each gene found, look within gene boundary for other read
                bool found = false;
                for (std::vector<GTFGene>::iterator it = results.begin(); it != results.end(); ++it) {
                    if (gtf->GetGene(it->GeneID()).CheckBoundary(align1->rname, align1->pos)) {
                        intragene_unannotated_splices.push_back(AlignmentPair(align0, align1, flag, true, is_backspliced)); 
                        found = true; 
                        break;
                    }
                }   

                if (found) {
                    continue;
                } else {

                    flag |= 1 << ALIGNED_SAME_CHR;
                    intrachromosomal_splices.push_back(AlignmentPair(align0, align1, flag, true, is_backspliced));
                    continue;    
                }        
            }         	
        }
    }

    //Now we go through each of the three sets, prioritizing the cis-gene model, as before
    if (intragene_unannotated_splices.size() > 0) {
        
        /*        
        //If this set of splices passes the quality filter
         if (ProcessSplices(intragene_unannotated_splices, minDiff)) {
        
             //We ignore this for now
             if (intragene_unannotated_splices[0].is_backspliced) {
                 gtf->IntrageneCircularSplice(intragene_unannotated_splices[0].align1->rname, 
                                              intragene_unannotated_splices[0].align1->pos, 
                                              intragene_unannotated_splices[0].align1->pos_end,
                                              intragene_unannotated_splices[0].align2->rname, 
                                              intragene_unannotated_splices[0].align2->pos, 
                                              intragene_unannotated_splices[0].align2->pos_end, 
                                              string(read->getId(), read->getIdLength()));   
             
             } else {   
             
                 gtf->IntrageneUnannotatedSplice(intragene_unannotated_splices[0].align1->rname, 
                                                 intragene_unannotated_splices[0].align1->pos, 
                                                 intragene_unannotated_splices[0].align1->pos_end,
                                                 intragene_unannotated_splices[0].align2->rname, 
                                                 intragene_unannotated_splices[0].align2->pos, 
                                                 intragene_unannotated_splices[0].align2->pos_end, 
                                                 string(read->getId(), read->getIdLength()));                              
             }
         } 
         */
          
    } else if (intrachromosomal_splices.size() > 0) {
        
        for (std::vector<AlignmentPair>::iterator it = intrachromosomal_splices.begin(); it != intrachromosomal_splices.end(); ++it) {

            gtf->IntrachromosomalSplice(it->align1->rname, 
                                        it->align1->pos, 
                                        it->align1->pos_end,
                                        it->align2->rname, 
                                        it->align2->pos, 
                                        it->align2->pos_end, 
                                        string(read->getId(), read->getIdLength()),
                                        string(read->getData(), read->getDataLength()));        

        }
    
    } else if (interchromosomal_splices.size() > 0) {
    
        for (std::vector<AlignmentPair>::iterator it = interchromosomal_splices.begin(); it != interchromosomal_splices.end(); ++it) {

            gtf->InterchromosomalSplice(it->align1->rname, 
                                        it->align1->pos, 
                                        it->align1->pos_end,
                                        it->align2->rname, 
                                        it->align2->pos, it->align2->pos_end,
                                        string(read->getId(), read->getIdLength()),
                                        string(read->getData(), read->getDataLength()));        
        }
    }
    
    for (std::vector<Alignment*>::iterator it = alignments.begin(); it != alignments.end(); ++it) {
        delete (*it);
    }
    
}

bool AlignmentFilter::ProcessSplices(std::vector<AlignmentPair> &pairs, unsigned minDiff) {

     if (pairs.size() == 1) {
        return true;
            
    } else {
         
        //Sort the scores by score, using operator< in AlignmentPair class
        //In this case we want to sort in the 'correct' order, so we reverse them after sorting
        sort(pairs.begin(), pairs.end());
        reverse(pairs.begin(), pairs.end());
                       
        //Check to see if the best alignment exceeds the second best alignment by at least confDiff
        unsigned diff = pairs[0].score - pairs[1].score;              
        if (diff >= minDiff) {
            return true; 
        } else {
            return false;       
        }  
    }   
}

void AlignmentFilter::FindPartialMatches(PairedAlignmentResult *result, AlignmentPair &pair) {

    seed_map map0, mapRC0, map1, mapRC1;
            
    //Temp variables, but may be able to use them later
    unsigned location = InvalidGenomeLocation;
    Direction direction;
    int score;
    int mapq;
    
    specialAligner->setReadId(0);
    specialAligner->CharacterizeSeeds(read0, map0, mapRC0);

    specialAligner->setReadId(1); 
    specialAligner->CharacterizeSeeds(read1, map1, mapRC1);
    
    //Print these maps
    //PrintMaps(map0, mapRC0);
    //PrintMaps(map1, mapRC1); 
    
    //The goal here is to look for possible partial alignments that could occur between
    //reads that are close together
    
    //FIX THIS 
    unsigned min_size = 1;
    
    //Loop over the map and mapRC for read1, adding any locations that pass the cutoff
    //into a vector of possible locations
    vector<unsigned> locs0, locs1;
    for (seed_map::iterator it = map0.begin(); it != map0.end(); ++it) {
        if (it->second.size() >= min_size) {
            locs0.push_back(it->first + *(it->second.begin()));
        }
    }
    
    for (seed_map::iterator it = mapRC0.begin(); it != mapRC0.end(); ++it) {
        if (it->second.size() >= min_size) {
            locs0.push_back(it->first + (read0->getDataLength() - *(it->second.rbegin())));
        }
    }
    
    for (seed_map::iterator it = map1.begin(); it != map1.end(); ++it) {
        if (it->second.size() >= min_size) {
            locs1.push_back(it->first + *(it->second.begin()));
        }
    }
    
    for (seed_map::iterator it = mapRC1.begin(); it != mapRC1.end(); ++it) {
        if (it->second.size() >= min_size) {
            locs1.push_back(it->first + (read1->getDataLength() - *(it->second.rbegin())));
        }
    }
       
    //Now loop over the possible locations, finding valid pairs
    for (vector<unsigned>::iterator it0 = locs0.begin(); it0 != locs0.end(); ++it0) {
        for (vector<unsigned>::iterator it1 = locs1.begin(); it1 != locs1.end(); ++it1) {

            const Genome::Contig *piece0 = genome->getContigAtLocation(*it0);
            string chr0 = piece0->name;
            int pos0 = (*it0) - piece0->beginningLocation + 1; 
            
            const Genome::Contig *piece1 = genome->getContigAtLocation(*it1);
            string chr1 = piece1->name;
            int pos1 = (*it1) - piece1->beginningLocation + 1; 
            
            if (chr0.compare(chr1) != 0) {
                continue;
            }
            
            //Calculate the distance between them
            unsigned distance = abs(pos1-pos0);            
            if (distance < maxSpacing) {
                result->status[0] = MultipleHits;
                result->status[1] = MultipleHits;
                result->mapq[0] = 1;
                result->mapq[1] = 1;
                return;
            }
        }    
    }    
}

void AlignmentFilter::CheckNoRC(PairedAlignmentResult *result, std::vector<AlignmentPair> &pairs) {

    //Go through the list of no_RC pairs and see if there is a match with the same chromosome and 
    //fewer mismatches
    for (std::vector<AlignmentPair>::iterator it = pairs.begin(); it != pairs.end(); ++it) {
    
        //If the chromosomes are the same for the noRC alignments
        if (it->align1->rname.compare(it->align2->rname) == 0) {
        
            //If the score is better
            if ((int) it->score < result->score[0] + result->score[1]) {
            
                result->status[0] = MultipleHits;
                result->status[1] = MultipleHits;
                result->mapq[0] = 1;
                result->mapq[1] = 1;
            }
        }
    }

}

void AlignmentFilter::ProcessPairs(PairedAlignmentResult* result, std::vector<AlignmentPair> &pairs) {

    if (pairs.size() == 1) {
           
		//First we convert the transcriptome coordinates back to a genomic position, for sorting
		if (pairs[0].align1->isTranscriptome) {

			result->tlocation[0] = pairs[0].align1->location;
			GenomeLocation offset;
			genome->getLocationOfContig(pairs[0].align1->rname.c_str(), &offset);
			offset += pairs[0].align1->pos - 1;
			result->location[0] = offset;

		} else {
			result->tlocation[0] = 0;
			result->location[0] = pairs[0].align1->location;
		}

		if (pairs[0].align2->isTranscriptome) {

			result->tlocation[1] = pairs[0].align2->location;
			GenomeLocation offset;
			genome->getLocationOfContig(pairs[0].align2->rname.c_str(), &offset);
			offset += pairs[0].align2->pos - 1;
			result->location[1] = offset;

		} else {
			result->tlocation[1] = 0;
			result->location[1] = pairs[0].align2->location;
		}

        //Only if the pair are both genome do we use the MAPQ value
        if (!pairs[0].align1->isTranscriptome && !pairs[0].align2->isTranscriptome) {
          genome_mapq = pairs[0].align1->mapq;
        }

        //Unique high quality hit
        result->status[0] = SingleHit;
        result->direction[0] = pairs[0].align1->direction;
        result->score[0] = pairs[0].align1->score;
        result->mapq[0] = min(maxMAPQ, genome_mapq);
        //result->mapq[0] = 70;
        result->isTranscriptome[0] = pairs[0].align1->isTranscriptome;
        
        result->status[1] = SingleHit;
        result->direction[1] = pairs[0].align2->direction;
        result->score[1] = pairs[0].align2->score;
        result->mapq[1] = min(maxMAPQ, genome_mapq);
        //result->mapq[1] = 70;
        result->isTranscriptome[1] = pairs[0].align2->isTranscriptome;
        
    } else {
         
        //Sort the scores by score, using operator< in AlignmentPair class
        sort(pairs.begin(), pairs.end());

		//First we convert the transcriptome coordinates back to a genomic position, for sorting
		if (pairs[0].align1->isTranscriptome) {

			result->tlocation[0] = pairs[0].align1->location;
			GenomeLocation offset;
			genome->getLocationOfContig(pairs[0].align1->rname.c_str(), &offset);
			offset += pairs[0].align1->pos - 1;
			result->location[0] = offset;

		} else {
			result->tlocation[0] = 0;
			result->location[0] = pairs[0].align1->location;
		}

		if (pairs[0].align2->isTranscriptome) {

			result->tlocation[1] = pairs[0].align2->location;
			GenomeLocation offset;
			genome->getLocationOfContig(pairs[0].align2->rname.c_str(), &offset);
			offset += pairs[0].align2->pos - 1;
			result->location[1] = offset;

		} else {
			result->tlocation[1] = 0;
			result->location[1] = pairs[0].align2->location;
		}

        //Only if the pair are both genome do we use the MAPQ value
        if (!pairs[0].align1->isTranscriptome && !pairs[0].align2->isTranscriptome) {
          genome_mapq = pairs[0].align1->mapq;
        }

        result->direction[0] = pairs[0].align1->direction;
        result->score[0] = pairs[0].align1->score;
        result->isTranscriptome[0] = pairs[0].align1->isTranscriptome;
        
        result->direction[1] = pairs[0].align2->direction;
        result->score[1] = pairs[0].align2->score;
        result->isTranscriptome[1] = pairs[0].align2->isTranscriptome;
    
        //Check to see if the best alignment exceeds the second best alignment
        //by at least confDiff
        unsigned diff = pairs[1].score - pairs[0].score;
                      
        if (diff >= confDiff) {
        
            //Unique high quality hit
            result->mapq[0] = min(maxMAPQ, genome_mapq);
            result->mapq[1] = min(maxMAPQ, genome_mapq);
            //result->mapq[0] = 70;
            //result->mapq[1] = 70;
            result->status[0] = SingleHit;
            result->status[1] = SingleHit;   
        
        } else {

            //Multiple hits
            result->mapq[0] = 1;
            result->mapq[1] = 1;
            result->status[0] = MultipleHits;
            result->status[1] = MultipleHits;           
        }  
    }
}

void AlignmentFilter::PrintMaps(seed_map& map, seed_map& mapRC) {

    printf("READ\n");
    for (seed_map::iterator it = map.begin(); it != map.end(); ++it) {
    
        printf("Pos: %u\n", it->first);
    
        const Genome::Contig *piece0 = genome->getContigAtLocation(it->first);
        const char* chr0 = piece0->name;
        unsigned pos0 = it->first - piece0->beginningLocation + 1; 
        
        printf("Pos: %s %u %u\n", chr0, pos0, it->second.size());
        for (std::set<unsigned>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            printf("%u\n", *it2);
        }   
    }
    
    printf("READRC\n");
    for (seed_map::iterator it = mapRC.begin(); it != mapRC.end(); ++it) {
    
        printf("Pos: %u\n", it->first);
    
        const Genome::Contig *piece0 = genome->getContigAtLocation(it->first);
        const char* chr0 = piece0->name;
        unsigned pos0 = it->first - piece0->beginningLocation + 1; 
        
        printf("Pos: %s %u %u\n", chr0, pos0, it->second.size());
        for (std::set<unsigned>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            printf("%u\n", *it2);
        }   
    }
}


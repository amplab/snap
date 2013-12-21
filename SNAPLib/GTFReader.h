/*++

Module Name:

    GTFReader.h

Abstract:

    Handles reading in GTF and GFF3 format annotation files. 

Authors:

    Andrew Magis, June, 2013

Environment:

Revision History:

--*/

#pragma once

//System includes
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "Compat.h"
#include "IntervalTree.h"
#include "Genome.h"
#include "FASTA.h"

#define GTF_MAX_READ_SIZE 4096

//Definitions for the possibilities of paired-end alignments          
const int FIRST_NOT_ALIGNED         = 0; 
const int SECOND_NOT_ALIGNED        = 1;
const int NOT_REVERSE_COMPLIMENTED  = 2;
const int ALIGNED_SAME_GENE         = 3;
const int ALIGNED_SAME_CHR          = 4;
const int ALIGNED_DIFF_CHR          = 5;
const int UNANNOTATED               = 6;
const int CIRCULAR                  = 7;

//Definitions for features types
enum {UNASSIGNED, EXON, INTRON};

class GTFReader;

//Namespaces
using namespace std;

typedef std::map<string,string> read_map;

class ReadInterval {

    friend class GTFReader;
    friend class ReadIntervalPair;
    friend class ReadIntervalMap;
    
    public:
        ReadInterval() : start(0), end(0), is_spliced(false), consolidated(false) {};
        ReadInterval(string chr, unsigned start_, unsigned end_, string id, string seq, bool is_spliced_);
        ReadInterval(const ReadInterval &rhs);
        virtual ~ReadInterval();
        ReadInterval& operator=(const ReadInterval &rhs);
        bool operator<(const ReadInterval &rhs) const;
        
        std::string Chr() const { return chr; };
        unsigned Start() const { return start; };
        unsigned End() const { return end; };
        std::string GeneID() const;
        std::string GeneName() const;
        std::string GeneNameSpliced(unsigned intersection) const;
        unsigned Intersection(const ReadInterval *rhs, read_map &difference) const;
        unsigned Difference(const ReadInterval *rhs, read_map &difference) const;
       
        void WriteGTF(ofstream &outfile, unsigned intersection) const;
        void Write(ofstream &outfile, unsigned intersection) const;
        void Write(ofstream &outfile) const;
        void WriteReads(ofstream &outfile) const;
        void Print() const;
        void GetGeneInfo(GTFReader *gtf);
    
    protected:
   
        bool Filter() const;
        void UpdateMatePointers(ReadInterval *new_pointer);
    
        string chr;
        unsigned start;
        unsigned end;
        read_map ids;
        set<string> gene_ids;
        set<string> gene_names;
        bool is_spliced;
        bool consolidated;
        bool counted;
        
        //Pointer to mate(s)
        set<ReadInterval*> mate;

};

class ReadIntervalPair {

    friend class GTFReader;
    friend class ReadIntervalMap;
    
    public:
        ReadIntervalPair(ReadInterval *interval1, ReadInterval *interval2);
        ReadIntervalPair(const ReadIntervalPair &rhs);
        virtual ~ReadIntervalPair();
        ReadIntervalPair& operator=(const ReadIntervalPair &rhs);
        bool operator<(const ReadIntervalPair &rhs) const;
        unsigned Intersection() const { return (unsigned) intersection.size(); };
 
        void Print() const;       
        void Write(ofstream &outfile, ofstream &readfile) const;
        void WriteGTF(ofstream &outfile) const;

    protected:
        
        ReadInterval *interval1;
        ReadInterval *interval2;
        read_map intersection;
        bool consolidated;
        
};

typedef std::pair<ReadIntervalPair, ReadIntervalPair> interval_pair;
typedef std::vector<interval_pair> spliced_mate;

class ReadIntervalMap {

    public:
        ReadIntervalMap();
        ReadIntervalMap(const ReadIntervalMap &rhs);
        virtual ~ReadIntervalMap();
        ReadIntervalMap& operator=(const ReadIntervalMap &rhs);
        
        void AddInterval(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq0, string seq1, bool is_spliced);
        void Consolidate(GTFReader *gtf, unsigned buffer, bool filterPromiscuousGenes);
        void Intersect(const ReadIntervalMap &map, unsigned buffer, unsigned minCount, GTFReader *gtf);
        
        void WriteSplicedMatePairs(ofstream &logfile, ofstream &readfile);
        void WriteGTF(ofstream &outfile);
        void Clear();
        void Print() const; 
        bool Find(GTFReader *gtf, string gene0, string gene1) const;       
        bool Find(GTFReader *gtf, ReadInterval *interval1, ReadInterval *interval2, string gene0, string gene1) const;

    protected:

        unsigned ConsolidateReadIntervals(unsigned buffer);
        ExclusiveLock mutex;
        
        //Vector of all paired-end reads that are not between genes
        std::vector<Interval<ReadInterval*> > read_intervals;
        IntervalTree<ReadInterval*> read_tree;  
        std::vector<ReadIntervalPair> pairs;  
        
        //Final results for a given map that contains overlaps between paired and spliced reads
        spliced_mate spliced_mate_pairs;

};

typedef std::set<string> id_set;

class GTFFeature {

    friend class GTFReader;
    friend class GTFTranscript;
    
    public:
    
        GTFFeature(string line);
        GTFFeature(const GTFFeature &rhs);
        virtual ~GTFFeature();
        GTFFeature& operator=(const GTFFeature &rhs);
        bool operator<(const GTFFeature &rhs) const;
        
	string Chr() const { return chr; };
        unsigned Start() const { return start; };
        unsigned End() const { return end; };
        unsigned ReadCount() const { return read_count; };
        unsigned Type() const { return type; };
        string TranscriptName() const;
        string GeneName() const;
        unsigned Length() const;
        bool GetAttribute(string key, string &value) const;
        void Print() const;
        
    protected:
    
        void IncrementReadCount();
    
        string key;
        string chr;
        string source;
        string feature;
        unsigned start;
        unsigned end;
        string score;
        char strand;
        char frame;
        id_set transcript_ids;
        unsigned type;

        string gene_id;
        string transcript_id;
        string gene_name;
        string transcript_name;
        std::map<std::string, std::string> attributes;
        unsigned read_count;
        ExclusiveLock mutex;

};

typedef std::map<std::string, GTFFeature> feature_map;
typedef std::map<std::string, GTFFeature*> feature_pointer_map;
typedef std::vector<GTFFeature*> feature_list;
typedef std::pair<unsigned, GTFFeature*> junction;

class GTFTranscript {

    friend class GTFReader;
    friend class GTFGene;    

    public:
        
        GTFTranscript(string chr, string gene_id, string transcript_id, string gene_name, string transcript_name, unsigned start, unsigned end);
        GTFTranscript(const GTFTranscript &rhs);
        virtual ~GTFTranscript();
        GTFTranscript& operator=(const GTFTranscript &rhs);
        
        string Chr() const { return chr; };
        string TranscriptID() const { return transcript_id; };
        string GeneID() const { return gene_id; };
        unsigned ReadCount() const { return (unsigned) read_count; };
        unsigned GenomicPosition(unsigned transcript_pos, unsigned span) const;
        void Junctions(unsigned start, unsigned span, std::vector<junction> &junctions) const;
        void WriteFASTA(const Genome *genome, std::ofstream &outfile) const;
        unsigned SplicedLength() const;        
        unsigned NormalizedCount() const;

        void Print() const;
        void WriteReadCountID(ofstream &outfile) const;
        void WriteReadCountName(ofstream &outfile) const;
        
    protected:
    
        void UpdateBoundaries(unsigned start, unsigned end);
        void Process(feature_map &all_features, feature_pointer_map &gene_features);
        void IncrementReadCount(unsigned numPotentialTranscripts);
    
        unsigned start;
        unsigned end;
        string chr;
        string gene_id;
        string gene_name;
        string transcript_id;
        string transcript_name;
        
        feature_list features;
        feature_list exons;
        float read_count;
        ExclusiveLock mutex;

};

typedef std::map<std::string, GTFTranscript> transcript_map;

class GTFGene {

    friend class GTFReader;
    
    public:
    
        GTFGene(string chr, string gene_id, unsigned start, unsigned end, string gene_name_);
        GTFGene(const GTFGene &rhs);
        virtual ~GTFGene();
        GTFGene & operator=(const GTFGene &rhs);
        bool operator<(const GTFGene &rhs) const;
        
        string Chr() const { return chr; };
        string GeneID() const { return gene_id; };
        string GeneName() const { return gene_name; };
        unsigned ReadCount() const { return read_count; };
        bool CheckBoundary(string query_chr, unsigned query_pos, unsigned buffer=1000) const;
        void Print() const;
        
        void WriteReadCountID(ofstream &outfile) const;
        void WriteReadCountName(ofstream &outfile) const;
        void WriteJunctionCountID(ofstream &outfile) const;     
   
    protected:
  
        void Process(feature_map &all_features, transcript_map &all_transcripts);  
        void UpdateBoundaries(unsigned start, unsigned end);
        void IncrementReadCount();
      

        string chr;
        string gene_id;
        string gene_name;
        unsigned start;
        unsigned end;
        feature_pointer_map features;
        id_set transcript_ids;
        unsigned read_count;
        ExclusiveLock mutex;
                          
};        

typedef std::map<std::string, GTFGene> gene_map;

class GTFReader {

    public:

        GTFReader(const char* output = NULL);
        virtual ~GTFReader();
        
        void Load(string filename);
        const GTFTranscript& GetTranscript(string transcript_id) const;
        const GTFGene& GetGene(string gene_id) const;
        unsigned Size() const { return (unsigned) transcripts.size(); };
        
        void IntervalGenes(std::string chr, unsigned start, unsigned stop, std::vector<GTFGene> &results);
        void IntervalTranscripts(std::string chr, unsigned start, unsigned stop, std::vector<GTFTranscript> &results);
        void IntervalFeatures(std::string chr, unsigned start, unsigned stop, std::vector<GTFFeature*> &results);
        
        //Function for incrementing read counts for each gene and transcript
        void IncrementReadCount(string transcript_id0, unsigned transcript_start0, unsigned start0, unsigned length0, string transcript_id1, unsigned transcript_start1, unsigned start1, unsigned length1);
        void IncrementReadCount(string transcript_id0, unsigned transcript_start0, unsigned start0, unsigned length0);        


        //Functions for building the transcriptome file
        void BuildTranscriptome(const Genome *genome);
        
        void IntrageneUnannotatedPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq0, string seq1);
        void IntrageneUnannotatedSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq);
        
        void IntrageneCircularPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq0, string seq1);
        void IntrageneCircularSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq);
        
        void IntrachromosomalPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq0, string seq1);
        void IntrachromosomalSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq);
        
        void InterchromosomalPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq0, string seq1);
        void InterchromosomalSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq);

        void AnalyzeReadIntervals();
        void WriteReadCounts();
               
        void Test();
    
    protected:

        int Parse(string line);
        unsigned ConsolidateReadIntervals(unsigned buffer);
        
        std::string filename;
        std::string prefix;
        feature_map features;
        transcript_map transcripts;
        gene_map genes;
        
        //Interval trees
        std::vector<Interval<GTFGene*> > gene_intervals;
        IntervalTree<GTFGene*> gene_tree;
                
        std::vector<Interval<GTFTranscript*> > transcript_intervals;
        IntervalTree<GTFTranscript*> transcript_tree;
    
        std::vector<Interval<GTFFeature*> > feature_intervals;
        IntervalTree<GTFFeature*> feature_tree;    
        
        //Reads that map within a single gene (but unannotated)
        ReadIntervalMap intragene_unannotated_pairs;
        ReadIntervalMap intragene_unannotated_splices;
        
        //Reads that map within a single gene (circularized)
        ReadIntervalMap intragene_circular_pairs;
        ReadIntervalMap intragene_circular_splices;        
        
        //Intrachromosomal reads
        ReadIntervalMap intrachromosomal_pairs;
        ReadIntervalMap intrachromosomal_splices;  
        
        //Intrachromosomal reads
        ReadIntervalMap interchromosomal_pairs;
        ReadIntervalMap interchromosomal_splices;           
        
};

inline std::string ToString(const unsigned& arg)
{
  char buffer[65];
  sprintf(buffer, "%u", arg) ;
  return std::string(buffer);
}



#pragma once

//System includes
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#define GTF_MAX_READ_SIZE 4096

//Namespaces
using namespace std;

class GTFFeature {

    friend class GTFReader;
    friend class GTFTranscript;
    
    public:
    
        GTFFeature(string line);
        GTFFeature(const GTFFeature &rhs);
        virtual ~GTFFeature();
        GTFFeature& operator=(const GTFFeature &rhs);
        bool operator<(const GTFFeature &rhs) const;
        
        unsigned Length() const;
        void Print() const;
        
        
    protected:
    
        string chr;
        string source;
        string feature;
        unsigned start;
        unsigned end;
        string score;
        char strand;
        char frame;
        string gene_id;
        string transcript_id;

};

typedef std::vector<GTFFeature> feature_list;
typedef std::pair<unsigned, unsigned> junction;

class GTFTranscript {

    friend class GTFReader;
    
    public:
        
        GTFTranscript(string chr, string gene_id, string transcript_id);
        GTFTranscript(const GTFTranscript &rhs);
        virtual ~GTFTranscript();
        GTFTranscript& operator=(const GTFTranscript &rhs);
        
        string Chr() { return chr; };
        string TranscriptID() { return transcript_id; };
        unsigned GenomicPosition(unsigned transcript_pos);
        std::vector<junction> Junctions(unsigned start, unsigned span);
        
    protected:
    
        void Process();
    
        string chr;
        string gene_id;
        string transcript_id;
        feature_list features;
        feature_list exons;

};

typedef std::map<std::string, GTFTranscript> transcript_map;

class GTFReader {

    public:

        GTFReader();
        virtual ~GTFReader();
        
        int Load(string filename);
        GTFTranscript& GetTranscript(string transcript_id);
        
        unsigned Size() { return transcripts.size(); };
        void Test();
    
    protected:

        int Parse(string line);
        string filename;
        transcript_map transcripts;
                
};


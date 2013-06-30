//Source includes
#include "GTFReader.h"

//System includes
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "Compat.h"

GTFFeature::GTFFeature(string line) {

    char _chr[GTF_MAX_READ_SIZE];
    char _source[GTF_MAX_READ_SIZE];
    char _feature[GTF_MAX_READ_SIZE];
    unsigned _start;
    unsigned _end;
    char _score[GTF_MAX_READ_SIZE];
    char _strand;
    char _frame;
    char _gene_id[GTF_MAX_READ_SIZE];
    char _transcript_id[GTF_MAX_READ_SIZE];
    
    //Tried stringstream - it is horrible, not thread safe, incredibly slow in a multithreaded application
    sscanf(line.c_str(), "%s\t%s\t%s\t%u\t%u\t%s\t%c\t%c\t%*[^\"]\"%255[^\"]\"%*[^\"]\"%255[^\"]\"", _chr, _source, _feature, &_start, &_end, _score, &_strand, &_frame, &_gene_id, &_transcript_id);
        
    chr = _chr;
    source = _source;
    feature = _feature;
    start = _start;
    end = _end;
    score = _score;
    strand = _strand;
    frame = _frame;
    gene_id = _gene_id;
    transcript_id = _transcript_id;
    
}

unsigned GTFFeature::Length() const {
    return (end-start)+1;
}

void GTFFeature::Print() const {
    printf("%s\t%s\t%s\t%u\t%u\t%s\t%c\t%c\t%s\t%s\n", chr.c_str(), source.c_str(), feature.c_str(), start, end, score.c_str(), strand, frame, gene_id.c_str(), transcript_id.c_str());
}

GTFFeature::GTFFeature(const GTFFeature& rhs) 
    :   chr(rhs.chr), source(rhs.source), feature(rhs.feature), start(rhs.start), end(rhs.end), score(rhs.score),
        strand(rhs.strand), frame(rhs.frame), gene_id(rhs.gene_id), transcript_id(rhs.transcript_id)
{}
  
GTFFeature::~GTFFeature() {}

GTFFeature& GTFFeature::operator=(const GTFFeature& rhs) {

    if (this != &rhs) {
        chr = rhs.chr;
        source = rhs.source;
        feature = rhs.feature;
        start = rhs.start;
        end = rhs.end;
        score = rhs.score;
        strand = rhs.strand;
        frame = rhs.frame;
        gene_id = rhs.gene_id;
        transcript_id = rhs.transcript_id;
    }
    return *this;
}

bool GTFFeature::operator<(const GTFFeature& rhs) const {
    return start < rhs.start;
}

GTFTranscript::GTFTranscript(string _chr, string _gene_id, string _transcript_id) 
    : chr(_chr), gene_id(_gene_id), transcript_id(_transcript_id) 
{}

GTFTranscript::GTFTranscript(const GTFTranscript& rhs) 
    : chr(rhs.chr), gene_id(rhs.gene_id), transcript_id(rhs.transcript_id), features(rhs.features)
{}
  
GTFTranscript::~GTFTranscript() {}

GTFTranscript& GTFTranscript::operator=(const GTFTranscript& rhs) {

    if (this != &rhs) {
        chr = rhs.chr;
        gene_id = rhs.gene_id;
        transcript_id = rhs.transcript_id;
        features = rhs.features;   
    }
    return *this;
}

void GTFTranscript::Process() {

    //Copy all exons into separate vector
    for (feature_list::iterator it = features.begin(); it != features.end(); ++it) {
        if (it->feature.compare("exon") == 0) {
            exons.push_back(*it);
        }
    }

    //Sort the exons by start and end
    sort(exons.begin(), exons.end());
    
}

unsigned GTFTranscript::GenomicPosition(unsigned transcript_pos) {
       
    //This assumes 1-offset transcript pos, and returns 1-offset genomic position
    //Converts transcript coordinates into genomic coordinates
    for (feature_list::iterator it = exons.begin(); it != exons.end(); ++it) {
        
        //If transcript_pos is less than or equal to this
        if (transcript_pos > it->Length()) {
            transcript_pos -= it->Length();
        } else {
            return it->start + transcript_pos - 1;
        }
    }
    printf("Warning, transcript_pos exceeds transcript length\n");
    return 0;
}

std::vector<junction> GTFTranscript::Junctions(unsigned transcript_pos, unsigned span) {

    //Container for splice junctions
    std::vector<junction> junctions;
    
    //Go through each feature of this transcript until we find the start
    unsigned current_pos = 0;
    feature_list::iterator current = exons.begin();
    for (feature_list::iterator next = ++(exons.begin()); next != exons.end(); ++next) {

        //Get the end of this exon
        current_pos += current->Length();

        //If transcript_pos is less than or equal to current, the start
        //lies within this feature
        if (transcript_pos <= current_pos) {

            //Check to see if the span exceeds the current feature. If so, 
            //simply return the current list of junctions
            if (current_pos-transcript_pos+1 >= span) {
                return junctions;
                
            //Otherwise, add the junction to the list of junctions
            } else {
            
                junctions.push_back(junction(current_pos+1, next->start-current->end-1));
                span -= current_pos-transcript_pos+1;
            
            }
        }
        current = next;            
    }
    return junctions;

}

//Constructor
GTFReader::GTFReader() {}

//Destructor
GTFReader::~GTFReader() {}

int GTFReader::Load(string _filename) {

    //Save the input filename
    filename = _filename;
    
    printf("Loading genome annotation from file... ");
    fflush(stdout);
    _int64 loadStart = timeInMillis();
    
    //Open input file
    std::ifstream infile(filename.c_str(), std::ios::in);
    if (!infile.is_open()) { 
        printf("Cannot open input file '%s'\n", filename.c_str());
	    exit(1);
    }
    
    string line;
    unsigned num_lines = 0;
    std::getline(infile, line, '\n');
    while (!infile.eof()){
  
        //Process this line
        Parse(line);
        num_lines++;
          
        std::getline(infile, line, '\n');
    }
    infile.close();
    
    //Sort each transcript
    for (transcript_map::iterator it = transcripts.begin(); it != transcripts.end(); ++it) {
        (*it).second.Process();
    }
    
    _int64 loadTime = timeInMillis() - loadStart;
    printf("%llds.  %u features,  %u transcripts\n", loadTime / 1000, num_lines, transcripts.size());
    
}

int GTFReader::Parse(string line) {

    // Create a new GTFFeature from this line
    GTFFeature feature(line);
    
    //Try to find this transcript in the transcript_map
    transcript_map::iterator pos = transcripts.find(feature.transcript_id);
        
    //If this sequence is not found, create a new vector to store this sequence (and others like it)
    if ((pos == transcripts.end())) {
        GTFTranscript transcript(feature.chr, feature.gene_id, feature.transcript_id);
        transcript.features.push_back(feature);
        transcripts.insert(transcript_map::value_type(feature.transcript_id, transcript));
        
    //Otherwise, add this feature to the transcript
    } else {
        pos->second.features.push_back(feature);
    }
    return 0;
    
}

GTFTranscript& GTFReader::GetTranscript(string transcript_id) {
 
    transcript_map::iterator pos = transcripts.find(transcript_id);
    if (pos == transcripts.end()) {
        //raise exception
        printf("No transcript %s\n", transcript_id.c_str());
        exit(1);
    }
    return pos->second;

}

void GTFReader::Test() {

    GTFTranscript &transcript = GetTranscript("ENST00000489673");

    std::vector<junction> junctions = transcript.Junctions(1, 200);
    
    for (std::vector<junction>::iterator it = junctions.begin(); it != junctions.end(); ++it) {
        printf("[%d %d]\n", it->first, it->second);
    }


}

/*
int main(int argc, char *argv[]) {

    const char *filename = NULL;

    //Process the command-line arguments
    int c;
    while ((c=getopt(argc, argv, "f:")) != EOF) {
    
        switch(c) {
            case 'f':
                filename = optarg;
                break;
            case ':':
                cerr << "Invalid option " << optarg << endl;
                exit(1);
                break; 	
	    }
    }
    
    struct timeval t1;
    gettimeofday(&t1, NULL);
    
    //Create the gtfReader
    GTFReader gtf;
    gtf.Load(filename);
    gtf.Test();
    
    struct timeval t2;
    gettimeofday(&t2, NULL);
    
    double elapsedTime;
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;
    fprintf(stderr, "Time elapsed: %.3lf seconds\n", elapsedTime / 1000.0);
    return 0;
}
*/


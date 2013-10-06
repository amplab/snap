/*++

Module Name:

    ContaminationFilter.cpp

Abstract:

    Handles counting of contaminant reads and printing output

Authors:

    Andrew Magis, June, 2013

Environment:

Revision History:

--*/

#include <vector>
#include "ContaminationFilter.h"

ContaminationFilter::ContaminationFilter(const Genome *contamination_, const char* output)
  : contamination(contamination_)
{
    prefix = "default";
    if (NULL != output) {
        prefix = output;
        size_t pos;
        if ((pos = prefix.rfind('.')) != std::string::npos) {
            prefix = prefix.substr(0, pos);
        }
    }
    pthread_mutex_init(&mutex, NULL);
}

ContaminationFilter::~ContaminationFilter() 
{
  pthread_mutex_destroy(&mutex);
}

int ContaminationFilter::AddAlignment(unsigned location, Direction direction, int score, int mapq, bool isTranscriptome, bool isMate0) {

    //Get the position and rname for this alignment
    string rname = "*";
    unsigned pos = 0;
    
    //If this is, in fact, aligned to something
    if (location != 0xFFFFFFFF) {
    
      const Genome::Piece *piece = contamination->getPieceAtLocation(location);
      rname = piece->name;
      pos = location - piece->beginningOffset + 1;
    
    } 
    
    //If the genomic location is valid
    if (pos != 0) {

       //Lock the mutex
       pthread_mutex_lock(&mutex);

       //Try to find this rname in the counts map
       contamination_count::iterator pos = counts.find(rname);
       
       //If it is not found, then add it
       if (pos == counts.end()) {
         counts.insert(contamination_count::value_type(rname, Contaminant(rname, 1)));

       //If it is found, then increment the count
       } else {
         pos->second.Increment();
       }

       //Unlock the mutex
       pthread_mutex_unlock(&mutex);

    }

}

void ContaminationFilter::Write() {

    ofstream contaminant_file;
    contaminant_file.open((prefix+".contaminants.txt").c_str());

    std::vector<Contaminant> temp;
    for (contamination_count::const_iterator it = counts.begin(); it != counts.end(); ++it) {
        temp.push_back(it->second);
    }
    sort(temp.rbegin(), temp.rend());

    for (std::vector<Contaminant>::const_iterator it = temp.begin(); it != temp.end(); ++it) {
      contaminant_file << it->Rname() << '\t' << it->Count() << endl;
    }
    contaminant_file.close();
}

void ContaminationFilter::Print() {

  for (contamination_count::const_iterator it = counts.begin(); it != counts.end(); ++it) {

    printf("%s: %u\n", it->second.Rname().c_str(), it->second.Count());


  } 




}

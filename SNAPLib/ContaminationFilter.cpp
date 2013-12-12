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
    InitializeExclusiveLock(&lock);
}

ContaminationFilter::~ContaminationFilter() 
{
    DestroyExclusiveLock(&lock);
}

int ContaminationFilter::AddAlignment(unsigned location, string header, string seq) {

    //Get the position and rname for this alignment
    string rname = "*";
    unsigned pos = 0;
    
    //If this is, in fact, aligned to something
    if (location != 0xFFFFFFFF) {
    
      const Genome::Contig *piece = contamination->getContigAtLocation(location);
      rname = piece->name;
      pos = location - piece->beginningOffset + 1;
    
    } 
    
    //If the genomic location is valid
    if (pos != 0) {

       //Lock the mutex
       AcquireExclusiveLock(&lock);

       //Get the existing element or insert a new one
       contamination_count::iterator pos = counts.insert(counts.begin(), contamination_count::value_type(rname, Contaminant(rname)));
       
       //Add the reads to the list
       pos->second.Add(header, seq);

       //Unlock the mutex
       ReleaseExclusiveLock(&lock);

    }

}

void ContaminationFilter::Write() {

    printf("Writing contaminants... ");
    fflush(stdout);
    _int64 loadStart = timeInMillis();

    ofstream contaminant_file, read_file;
    contaminant_file.open((prefix+".contaminants.txt").c_str());
    read_file.open((prefix+".contaminants.reads.fa").c_str());

    //Sort the contaminants by the number of reads there are
    std::vector<Contaminant> temp;
    for (contamination_count::const_iterator it = counts.begin(); it != counts.end(); ++it) {
        temp.push_back(it->second);
    }
    sort(temp.rbegin(), temp.rend());

    //Go through each one and output the reads
    for (std::vector<Contaminant>::const_iterator it = temp.begin(); it != temp.end(); ++it) {
      contaminant_file << it->Rname() << '\t' << it->Count() << endl;
      it->WriteReads(read_file);
    }
    contaminant_file.close();
    read_file.close();

    _int64 loadTime = timeInMillis() - loadStart;
    printf("%llds.\n", loadTime / 1000);
    fflush(stdout);

}

void ContaminationFilter::Print() {

  for (contamination_count::const_iterator it = counts.begin(); it != counts.end(); ++it) {

    printf("%s: %u\n", it->second.Rname().c_str(), it->second.Count());


  } 




}

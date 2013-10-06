/*++

Module Name:

    ContaminationFilter.h

Abstract:

    Handles counting of contaminant reads and printing output

Authors:

    Andrew Magis, June, 2013

Environment:

Revision History:

--*/

#pragma once

//System includes
#include <string>
#include <stdio.h>
#include <map>
#include <fstream>
#include <pthread.h>

#include "PairedEndAligner.h"
#include "Genome.h"

using namespace std;

class Contaminant {

    public:

      Contaminant(std::string rname_, unsigned count_=0) : rname(rname_), count(count_) {};
      Contaminant(const Contaminant &rhs) : rname(rhs.rname), count(rhs.count) {};
      Contaminant& operator=(const Contaminant &rhs) { rname = rhs.rname; count = rhs.count; };
      bool operator<(const Contaminant &rhs) const { return count < rhs.count; };
      virtual ~Contaminant() {};
      
      void Increment() { count++; };
      unsigned Count() const { return count; };
      std::string Rname() const { return rname; };

    private:

      std::string rname;
      unsigned count;
};

typedef std::map<std::string, Contaminant> contamination_count;

class ContaminationFilter {
  
    public:
    
        //Constructors/destructor
        ContaminationFilter(const Genome* contamination, const char* output);
        virtual ~ContaminationFilter();  
        
        //Functions
        int AddAlignment(unsigned location, Direction direction, int score, int mapq, bool isTranscriptome, bool isMate0); 
        void Print();     
        void Write();

    protected:
        
        const Genome* contamination;
        contamination_count counts;
        pthread_mutex_t mutex;
        string prefix;
};



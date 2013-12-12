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

#include "stdafx.h"
#include "Compat.h"
#include "PairedEndAligner.h"
#include "Genome.h"

using namespace std;

typedef std::map<string, string> read_map;

class Contaminant {

    public:

      Contaminant(std::string rname_) : rname(rname_) {};
      Contaminant(const Contaminant &rhs) : rname(rhs.rname), reads(rhs.reads) {};
      Contaminant& operator=(const Contaminant &rhs) { rname = rhs.rname; reads = rhs.reads; };
      bool operator<(const Contaminant &rhs) const { return Count() < rhs.Count(); };
      virtual ~Contaminant() {};
      
      void Add(string header, string seq) { reads.insert(read_map::value_type(header, seq)); };
      unsigned Count() const { return reads.size(); };
      std::string Rname() const { return rname; };

      void WriteReads(ofstream &outfile) const {
        for (read_map::const_iterator it = reads.begin(); it != reads.end(); ++it) {
            outfile << ">"+it->first+"|"+rname << endl << it->second << endl;
        }
      };

    private:

      std::string rname;
      read_map reads;
};

typedef std::map<std::string, Contaminant> contamination_count;

class ContaminationFilter {
  
    public:
    
        //Constructors/destructor
        ContaminationFilter(const Genome* contamination, const char* output);
        virtual ~ContaminationFilter();  
        
        //Functions
        int AddAlignment(unsigned location, string header, string seq); 
        void Print();     
        void Write();

    protected:
        
        const Genome* contamination;
        contamination_count counts;
        ExclusiveLock lock;
        string prefix;
};



/*++

Module Name:

    GTFReader.cpp

Abstract:

    Handles reading in GTF and GFF3 format annotation files. 

Authors:

    Andrew Magis, June, 2013

Environment:

Revision History:

--*/

#include "stdafx.h"

//System includes
#include <stdio.h>
#include <time.h>
#include <string.h>

//Source includes
#include "GTFReader.h"
#include "Compat.h"
#include "exit.h"

using std::min;
using std::max;

bool IntervalNodeSort(const Interval<ReadInterval*> &interval0, const Interval<ReadInterval*> &interval1) { 
    return interval0.value->Start() < interval1.value->Start(); 
}

bool FeatureListSort(const GTFFeature* feature0, const GTFFeature* feature1) {
    return feature0->Start() < feature1->Start();
}

bool SplicedMateSort(const interval_pair &pair0, const interval_pair &pair1) {
    return ((pair0.first.Intersection() + pair0.second.Intersection()) >
            (pair1.first.Intersection() + pair1.second.Intersection()));

}

ReadInterval::ReadInterval(string chr_, unsigned start_, unsigned end_, string id, string seq, bool is_spliced_) 
    : chr(chr_), start(start_), end(end_), is_spliced(is_spliced_), consolidated(false), counted(false)
{
    ids.insert(read_map::value_type(id, seq));
} 
    
ReadInterval::ReadInterval(const ReadInterval &rhs) 
    : chr(rhs.chr), start(rhs.start), end(rhs.end), ids(rhs.ids), is_spliced(rhs.is_spliced), 
      gene_ids(rhs.gene_ids), gene_names(rhs.gene_names), consolidated(false), counted(false), mate(rhs.mate)
{}

ReadInterval::~ReadInterval() 
{}

ReadInterval& ReadInterval::operator=(const ReadInterval &rhs) {
    if (this != &rhs) {
        chr = rhs.chr;
        start = rhs.start;
        end = rhs.end;
        ids = rhs.ids;
        gene_ids = rhs.gene_ids;
        gene_names = rhs.gene_names;
        is_spliced = rhs.is_spliced;
        consolidated = rhs.consolidated;
        counted = rhs.counted;
        mate = rhs.mate;
    }
    return *this;    
}

void ReadInterval::UpdateMatePointers(ReadInterval *new_pointer) {
    
    //Go through each mate
    for (set<ReadInterval*>::iterator it = mate.begin(); it != mate.end(); ++it) {
    
        //Delete the keys that point to this object
        (*it)->mate.erase(this); 
        
        //Replace it with the new pointer
        (*it)->mate.insert(new_pointer);
    
    }
}

unsigned ReadInterval::Intersection(const ReadInterval *rhs, read_map &intersection) const {

    //Return the size of the intersection between the two sets of IDs
    for (read_map::const_iterator it = ids.begin(); it != ids.end(); ++it) {
              
        read_map::const_iterator pos = rhs->ids.find(it->first);
        if (pos != rhs->ids.end()) {
            intersection.insert(*it);
        }        
    }
    return (unsigned) intersection.size();
}

unsigned ReadInterval::Difference(const ReadInterval *rhs, read_map &difference) const {

    //Return the size of the intersection between the two sets of IDs
    for (read_map::const_iterator it = ids.begin(); it != ids.end(); ++it) {
        
        read_map::const_iterator pos = rhs->ids.find(it->first);
        if (pos == rhs->ids.end()) {
            difference.insert(*it);
        }        
    }
    return (unsigned) difference.size();
}

bool ReadInterval::operator<(const ReadInterval& rhs) const {
    return (start <= rhs.start);
}

std::string ReadInterval::GeneID() const {

    std::set<string>::iterator it = gene_ids.begin();
    if (it == gene_ids.end()) {
        return "NoGene";
    }
    string gene_id = (*it);
    for (++it; it != gene_ids.end(); ++it) {
        gene_id += ',' + (*it);
    }
    return gene_id;
}

std::string ReadInterval::GeneName() const {

    std::set<string>::iterator it = gene_names.begin();
    if (it == gene_names.end()) {
        return GeneID();
    }
    string gene_name = (*it);
    for (++it; it != gene_names.end(); ++it) {
        gene_name += ',' + (*it);
    }
    return gene_name;
}

std::string ReadInterval::GeneNameSpliced(unsigned intersection) const {

    //Get gene_name, if present or gene_id if not
    string gene_name = GeneName();
    if (is_spliced) {
        gene_name += ",S";
    } else {
        gene_name += ",P";
    }
    return gene_name + "," + ToString(intersection);
}

void ReadInterval::GetGeneInfo(GTFReader *gtf) {

    //Get gene info for this interval
    std::vector<GTFGene> results;
    gtf->IntervalGenes(chr, start, end, results);
    
    for (std::vector<GTFGene>::iterator it = results.begin(); it != results.end(); ++it) {

        gene_ids.insert(it->GeneID());
        if (it->GeneName().size() > 0) {
            gene_names.insert(it->GeneName());
        }
    }

}

bool ReadInterval::Filter() const {

   if (chr.find("MT") != std::string::npos) {
      return true;
   } 

   if (chr.find("GL") != std::string::npos) {
      return true;
   }

   for (std::set<string>::iterator it = gene_names.begin(); it != gene_names.end(); ++it) {
      if (it->find("HLA-") != std::string::npos) {
         return true;
      }
   }
   return false;

}

void ReadInterval::Write(ofstream &outfile, unsigned intersection) const {
    outfile << chr << ':' << start << '-' << end << '\t' << GeneID() << '\t' << GeneNameSpliced(intersection) << endl;
}

void ReadInterval::Write(ofstream &outfile) const {
    outfile << chr << ':' << start << '-' << end << '\t';
}

void ReadInterval::WriteReads(ofstream &outfile) const {
 
    string type = "splice";
    if (!is_spliced) {
        type = "mate";
    }

    //Write reads to read file
    for (read_map::const_iterator it = ids.begin(); it != ids.end(); ++it) {
        //Only write out reads that have a finite length
        if (it->second.size() > 0) {
            outfile << ">" << it->first << "|" << GeneName() << "|" << type << endl;
            outfile << it->second << endl;
        }
    }
}

void ReadInterval::Print() const {
    cout << chr << ':' << start << '-' << end << '\t' << GeneID() << '\t' << GeneName() << '\t' << ids.size() << '\t' << is_spliced;
    printf("\t%p [", this);
    for (std::set<ReadInterval*>::iterator it = mate.begin(); it != mate.end(); ++it) {
      printf("\t%p",(*it));
    }
    printf("]\n");
}

void ReadInterval::WriteGTF(ofstream &outfile, unsigned intersection) const {
    outfile << chr << '\t' << "snapr" << '\t' << "interval" << '\t' << start << '\t' << end << '\t' << '.' << '\t' << '.' << '\t' << '.' << '\t';
    outfile << "gene_id \"" << GeneID() << "\"; transcript_id \"" << GeneNameSpliced(intersection) << "\"; gene_name \"" << GeneName() << "\";" << endl;
}

ReadIntervalPair::ReadIntervalPair(ReadInterval *interval1_, ReadInterval *interval2_) 
    : interval1(interval1_), interval2(interval2_)
{
    interval1->Intersection(interval2, intersection);
} 
    
ReadIntervalPair::ReadIntervalPair(const ReadIntervalPair &rhs) 
    : interval1(rhs.interval1), interval2(rhs.interval2), intersection(rhs.intersection)
{}

ReadIntervalPair::~ReadIntervalPair() 
{}

ReadIntervalPair& ReadIntervalPair::operator=(const ReadIntervalPair &rhs) {
    if (this != &rhs) {
        interval1 = rhs.interval1;
        interval2 = rhs.interval2;
        intersection = rhs.intersection;
    }
    return *this;    
}

bool ReadIntervalPair::operator<(const ReadIntervalPair& rhs) const {
    return intersection.size() > rhs.intersection.size();
}

void ReadIntervalPair::Write(ofstream &outfile, ofstream &readfile) const {

    outfile << intersection.size() << '\t';
    interval1->Write(outfile);
    outfile << interval1->GeneID() << '\t' << interval1->GeneName() << '\t';
    interval2->Write(outfile);
    outfile << interval2->GeneID() << '\t' << interval2->GeneName();

    string type1 = "splice";
    if (!interval1->is_spliced) {
        type1 = "mate";
    }
    string type2 = "splice";
    if (!interval2->is_spliced) {
        type2 = "mate";
    }

    read_map reads;
    for (read_map::const_iterator it = interval1->ids.begin(); it != interval1->ids.end(); ++it) {
        //Add all reads to the read_map
        if (it->second.size() > 0) {
            reads.insert(read_map::value_type(it->first+'|'+type1+'|'+interval1->GeneName(), it->second));
        }
    }    
    for (read_map::const_iterator it = interval2->ids.begin(); it != interval2->ids.end(); ++it) {
        //Add all reads to the read_map
        if (it->second.size() > 0) {
            reads.insert(read_map::value_type(it->first+'|'+type2+'|'+interval2->GeneName(), it->second));
        }
    }

    //Write the sorted list to the output file
    for (read_map::const_iterator it = reads.begin(); it != reads.end(); ++it) {
        readfile << ">" << it->first << endl << it->second << endl;
    }    


    //interval1->WriteReads(readfile);
    //interval2->WriteReads(readfile);   
}

void ReadIntervalPair::Print() const {
    interval1->Print();
    interval2->Print();
}

void ReadIntervalPair::WriteGTF(ofstream &outfile) const {
    interval1->WriteGTF(outfile, (unsigned) intersection.size());
    interval2->WriteGTF(outfile, (unsigned) intersection.size());
}

ReadIntervalMap::ReadIntervalMap() {   
    InitializeExclusiveLock(&mutex);
}
   
ReadIntervalMap::ReadIntervalMap(const ReadIntervalMap &rhs) {
    printf("ReadIntervalMap: copy constructor not implemented\n");
    exit(1);
}

ReadIntervalMap::~ReadIntervalMap() {
    DestroyExclusiveLock(&mutex);

    for (std::vector<Interval<ReadInterval*> >::iterator it = read_intervals.begin(); it != read_intervals.end(); ++it) {
        delete it->value;
    }
    read_intervals.clear();
}

void ReadIntervalMap::Clear() {
    for (std::vector<Interval<ReadInterval*> >::iterator it = read_intervals.begin(); it != read_intervals.end(); ++it) {
        delete it->value;
    }
    read_intervals.clear();
}

ReadIntervalMap& ReadIntervalMap::operator=(const ReadIntervalMap &rhs) {
    printf("ReadIntervalMap: assignment operator not implemented\n");
    exit(1);
}

void ReadIntervalMap::AddInterval(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq0, string seq1, bool is_spliced) {

    //Get the mutex
    AcquireExclusiveLock(&mutex);
    
    ReadInterval *mate0 = new ReadInterval(chr0, start0, end0, id, seq0, is_spliced);
    ReadInterval *mate1 = new ReadInterval(chr1, start1, end1, id, seq1, is_spliced);

    mate0->mate.insert(mate1);
    mate1->mate.insert(mate0);
    
    //Add this to the vector of links
    read_intervals.push_back(Interval<ReadInterval*>(start0, end0, mate0));
    read_intervals.push_back(Interval<ReadInterval*>(start1, end1, mate1));
    
    //Unlock it
    ReleaseExclusiveLock(&mutex);
    
}

void ReadIntervalMap::Print() const {

    printf("Begin Map\n");
    for (std::vector<Interval<ReadInterval*> >::const_iterator it = read_intervals.begin(); it != read_intervals.end(); ++it) {
      it->value->Print();
    }
    printf("End Map\n");

}

void ReadIntervalMap::Consolidate(GTFReader *gtf, unsigned buffer, bool filterPromiscuousGenes=true) {

    //printf("Building interval tree of read pairs\n");
    //Print();

    //Get size of current interval set
    unsigned initial_size = (unsigned) read_intervals.size();
       
    //Repeatedly consolidate existing reads until it no longer changes
    do {
        initial_size = (unsigned) read_intervals.size();
        ConsolidateReadIntervals(buffer);
        //printf("Initial Size: %u Current Size: %u\n", initial_size, read_intervals.size());
    } while (initial_size > read_intervals.size());

    //Print();  

    std::vector<Interval<ReadInterval*> > temp_intervals;

    for (std::vector<Interval<ReadInterval*> >::iterator it = read_intervals.begin(); it != read_intervals.end(); ++it) {

       //Get inveral information
       it->value->GetGeneInfo(gtf);
        if (filterPromiscuousGenes && it->value->Filter()) {
           continue;
        }
       
        temp_intervals.push_back(Interval<ReadInterval*>(it->value->start, it->value->end, it->value)); 

    }

    //Build new tree with filtered intervals
    read_intervals = temp_intervals;
    //sort(read_intervals.begin(), read_intervals.end(), IntervalNodeSort);
    read_tree = IntervalTree<ReadInterval*>(read_intervals);

    //Now we build all ReadIntervalPairs
    for (std::vector<Interval<ReadInterval*> >::iterator it = read_intervals.begin(); it != read_intervals.end(); ++it) {
        
        //Get interval info
        it->value->GetGeneInfo(gtf);
     
        if (filterPromiscuousGenes && it->value->Filter()) {
           continue;
        }
                  
        //For each mate in this result
        for (std::set<ReadInterval*>::iterator mate_it = it->value->mate.begin(); mate_it != it->value->mate.end(); ++mate_it) {    
        
            //Now doing this in Intersect because info wasn't gettting though, and I never figured out why
            (*mate_it)->GetGeneInfo(gtf);

            if (filterPromiscuousGenes && (*mate_it)->Filter()) {
              continue;
            }

            //it->value->Print();
            //(*mate_it)->Print();

            //Create a new pair
            pairs.push_back(ReadIntervalPair(it->value, *mate_it)); 

            //Find(gtf, it->value, *mate_it, "EGFR", "PSPH");

            //Remove me from the linking, so that mate won't connect to me later
            (*mate_it)->mate.erase(it->value);
        }
        
    }
    
    //Sort the pairs by the number of mate pairs they have in common
    sort(pairs.begin(), pairs.end());
  
}

unsigned ReadIntervalMap::ConsolidateReadIntervals(unsigned buffer) {

    //Sort the read intervals, which avoids weird issue in IntervalTree
    //sort(read_intervals.begin(), read_intervals.end(), IntervalNodeSort);

    //Build the interval tree from the current set of intervals
    read_tree = IntervalTree<ReadInterval*>(read_intervals);
    
    //Create a new, temporary set of intervals
    std::vector<Interval<ReadInterval*> > temp_intervals;
    
    //Loop over all the read_intervals that we currently have
    for (std::vector<Interval<ReadInterval*> >::iterator it = read_intervals.begin(); it != read_intervals.end(); ++it) {

        //Check to see if start is 0, which means it has already been consolidated
        if (it->value->consolidated) {
            continue;
        }
    
        //Create a new interval the same as the current one
        ReadInterval *new_interval = new ReadInterval(*(it->value));

        //Update the mates to point at the new interval
        new_interval->UpdateMatePointers(new_interval);

        //Set myself as consolidated
        //it->value->consolidated = true;

        //Query the interval tree for this interval to find any overlaps
        std::vector<Interval<ReadInterval*> > temp_results;
        std::vector<ReadInterval*> chr_results;
        read_tree.findOverlapping(it->value->start-buffer, it->value->end+buffer, temp_results);
 
        //Filter by chromosome
        for (std::vector<Interval<ReadInterval*> >::iterator it2 = temp_results.begin(); it2 != temp_results.end(); ++it2) {
        
            //Make sure this has not already been consolidated
            if (it2->value->consolidated) {
                continue;
            }
        
            //Make sure the chromosomes match
            if (it2->value->chr.compare(it->value->chr) == 0) {
                chr_results.push_back(it2->value);
            }   
        }   

        //I have already been consolidated
        if (chr_results.size() == 0) {
            //printf("Warning, cannot find myself\n");
        } 

        for (std::vector<ReadInterval*>::iterator it3 = chr_results.begin(); it3 != chr_results.end(); ++it3) {
            
            new_interval->start = min(new_interval->start, (*it3)->start);
            new_interval->end = max(new_interval->end, (*it3)->end);
            std::copy((*it3)->ids.begin(), (*it3)->ids.end(), std::inserter(new_interval->ids, new_interval->ids.end()));
       
            //Update the boundaries of any intervals that pointed to this interval
            (*it3)->UpdateMatePointers(new_interval);
            
            //Update the mates
            for (set<ReadInterval*>::iterator it4 = (*it3)->mate.begin(); it4 != (*it3)->mate.end(); ++it4) {
                new_interval->mate.insert(*it4);
            }
            
            //Indicate this interval has been used
            (*it3)->consolidated = true;
            
        }
                
        //Now we have built new interval. Insert it in the new interval vector
        temp_intervals.push_back(Interval<ReadInterval*>(new_interval->start, new_interval->end, new_interval));
    
    } 
    
//     //Delete all the original intervals that were zeroed out
//     for (std::vector<Interval<ReadInterval*> >::iterator it = read_intervals.begin(); it != read_intervals.end(); ++it) {
//         if (it->value->consolidated) {
//             delete it->value;
//         }
//     }
    
    
    //Finally, assign set of new intervals to the old set, and return the size
   
    read_intervals = temp_intervals;
    return (unsigned) read_intervals.size();

}

bool ReadIntervalMap::Find(GTFReader *gtf, string gene0, string gene1) const {

    //Go through each pair
    for (std::vector<ReadIntervalPair>::const_iterator it = pairs.begin(); it != pairs.end(); ++it) {
        
        it->interval1->GetGeneInfo(gtf);
        it->interval2->GetGeneInfo(gtf);
  
        std::set<string>::const_iterator pos1, pos2;
        if (((pos1 = it->interval1->gene_names.find(gene0)) != it->interval1->gene_names.end()) &&
            ((pos2 = it->interval2->gene_names.find(gene1)) != it->interval2->gene_names.end())) {
 
            printf("Found IntervalPair %s and %s\n", gene0.c_str(), gene1.c_str());
            it->Print();
            printf("\n");     
            return true;

        }

        if (((pos1 = it->interval1->gene_names.find(gene1)) != it->interval1->gene_names.end()) &&
            ((pos2 = it->interval2->gene_names.find(gene0)) != it->interval2->gene_names.end())) {
      
              printf("Found IntervalPair %s and %s\n", gene1.c_str(), gene0.c_str());
              it->Print();
              printf("\n");  
              return true;
          
         }
    }

    return false;
}

bool ReadIntervalMap::Find(GTFReader *gtf, ReadInterval *interval1, ReadInterval *interval2, string gene0, string gene1) const {

       interval1->GetGeneInfo(gtf);
       interval2->GetGeneInfo(gtf);

       std::set<string>::const_iterator pos1, pos2;
       if (((pos1 = interval1->gene_names.find(gene0)) != interval1->gene_names.end()) &&
           ((pos2 = interval2->gene_names.find(gene1)) != interval2->gene_names.end())) {

            printf("SMALL Found IntervalPair %s and %s\n", gene0.c_str(), gene1.c_str());
            interval1->Print();
            interval2->Print();
            printf("\n");
            return true;

        }

        if (((pos1 = interval1->gene_names.find(gene1)) != interval1->gene_names.end()) &&
            ((pos2 = interval2->gene_names.find(gene0)) != interval2->gene_names.end())) {

              printf("SMALL Found IntervalPair %s and %s\n", gene1.c_str(), gene0.c_str());
              interval1->Print();
              interval2->Print();
              printf("\n");
              return true;

         }

        return false;

}


void ReadIntervalMap::Intersect(const ReadIntervalMap &rhs, unsigned buffer, unsigned minCount, GTFReader *gtf) {


    //Loop over all intervals in the other map
    for (std::vector<ReadIntervalPair>::const_iterator it = rhs.pairs.begin(); it != rhs.pairs.end(); ++it) {

        std::vector<Interval<ReadInterval*> > left_temp_results;
        std::vector<ReadInterval*> left_chr_results;
        read_tree.findOverlapping(it->interval1->start-buffer, it->interval1->end+buffer, left_temp_results);

        //Filter by chromosome
        for (std::vector<Interval<ReadInterval*> >::iterator it2 = left_temp_results.begin(); it2 != left_temp_results.end(); ++it2) {
        
            //Make sure the chromosomes match
            if (it2->value->chr.compare(it->interval1->chr) == 0) {
                left_chr_results.push_back(it2->value);
            }   
        }
 
        std::vector<Interval<ReadInterval*> > right_temp_results;
        std::vector<ReadInterval*> right_chr_results;
        read_tree.findOverlapping(it->interval2->start-buffer, it->interval2->end+buffer, right_temp_results);

        //Filter by chromosome
        for (std::vector<Interval<ReadInterval*> >::iterator it2 = right_temp_results.begin(); it2 != right_temp_results.end(); ++it2) {
  
            //Make sure the chromosomes match
            if (it2->value->chr.compare(it->interval2->chr) == 0) {
               right_chr_results.push_back(it2->value);
            }   
        }     

        //Now all we have to do is verify that the results are linked to each other
        for (std::vector<ReadInterval*>::iterator left = left_chr_results.begin(); left != left_chr_results.end(); ++left) {
            for (std::vector<ReadInterval*>::iterator right = right_chr_results.begin(); right != right_chr_results.end(); ++right) {
        
                //Find(gtf, *left, *right, "EGFR", "PSPH");

                set<ReadInterval*>::iterator pos = (*left)->mate.find(*right);

                //If the mate of one contains the other
                if ((pos != (*left)->mate.end())) {
 
                    ReadIntervalPair pair0 = ReadIntervalPair(it->interval1, it->interval2);
                    ReadIntervalPair pair1 = ReadIntervalPair(*left, *right);
                    
                    //Check that the size of each pair exceeds the count
                    if ((pair0.intersection.size() >= minCount) && (pair1.intersection.size() >= minCount)) {

                        it->interval1->GetGeneInfo(gtf);
                        it->interval2->GetGeneInfo(gtf);
                        (*left)->GetGeneInfo(gtf);
                        (*right)->GetGeneInfo(gtf);

                        //Create a new ReadIntervalPair with both of these overlapping sets
                        spliced_mate_pairs.push_back( std::make_pair(ReadIntervalPair(it->interval1, it->interval2), ReadIntervalPair(*left, *right)) );

                    }
               } 
            }
        }
    }
}

void ReadIntervalMap::WriteGTF(ofstream &outfile) {

    //Sort the pairs
    sort(spliced_mate_pairs.begin(), spliced_mate_pairs.end(), SplicedMateSort);

    //Write out each one to the output file
    for (spliced_mate::const_iterator it = spliced_mate_pairs.begin(); it != spliced_mate_pairs.end(); ++it) {
    
        it->first.WriteGTF(outfile);
        it->second.WriteGTF(outfile);
       
    }
}

void ReadIntervalMap::WriteSplicedMatePairs(ofstream &logfile, ofstream &readfile) {

    //Sort the pairs
    sort(spliced_mate_pairs.begin(), spliced_mate_pairs.end(), SplicedMateSort);

    //Write out each one to the output file
    for (spliced_mate::const_iterator it = spliced_mate_pairs.begin(); it != spliced_mate_pairs.end(); ++it) {
     
        logfile << "Mated" << '\t';  
        it->first.Write(logfile, readfile);
        logfile << endl;
        logfile << "Spliced" << '\t';
        it->second.Write(logfile, readfile);
        logfile << endl << endl;
        
    }
}

GTFFeature::GTFFeature(string line) 
    : read_count(0) 
{

    InitializeExclusiveLock(&mutex);

    char* line_c = (char*)line.c_str(); 
    char *pch;
    pch = strtok(line_c,"'\t'"); chr = pch;
    key = string(pch);
    pch = strtok(NULL,"'\t'"); source = pch;
    pch = strtok(NULL,"'\t'"); feature = pch;
    pch = strtok(NULL,"'\t'"); start = atoi(pch);
    key += string(pch);
    pch = strtok(NULL,"'\t'"); end = atoi(pch);
    key += string(pch);
    pch = strtok(NULL,"'\t'"); score = pch;
    pch = strtok(NULL,"'\t'"); strand = *pch;
    pch = strtok(NULL,"'\t'"); frame = *pch;
       
    while (pch != NULL) {
    
        char *temp_key, *temp_value;
        pch = strtok(NULL, " ="); temp_key = pch;
        pch = strtok(NULL, ";"); temp_value = pch;
        
        if (temp_key != NULL) {
            
            string key = temp_key;
            string value = temp_value;
            
            //Remove any quotes or spaces from each string
            value.erase(remove(value.begin(), value.end(), '\"' ), value.end());
            attributes.insert(std::map<string,string>::value_type(key, value));
        
        }
    }

    //Set the type
    if (feature.compare("exon") == 0) {
        type = EXON;
    }
    
    string value;
    
    //If gene_id is present, use this as gene_id
    if (GetAttribute("gene_id", value)) {
        gene_id = value;
        
    //Otherwise, use Parent
    } else if (GetAttribute("Parent", value)) {
        gene_id = value;
    
    } else {
        gene_id = "Unknown";
    } 
    
    //If transcript_id exists
    if (GetAttribute("transcript_id", value)) {
        transcript_id = value;
    } else {
        transcript_id = gene_id;
    }

    //Prepend gene_id onto the feature name to prevent overlapping genes from conflicting
    key = gene_id + key;

}

bool GTFFeature::GetAttribute(string key, string &value) const {

    //At the end we must define the gene_id and transcript_id
    std::map<string,string>::const_iterator pos = attributes.find(key);
    if (pos == attributes.end()) {
        return false;
    }
    
    //Return the attribute
    value = pos->second;
    return true;
}

string GTFFeature::TranscriptName() const {
    string value;
    //if transcript_name exists
    if (GetAttribute("transcript_name", value)) {
        return value;
    } else {
        return transcript_id;
    }
}

string GTFFeature::GeneName() const {
    string value;
    //if gene_name exists
    if (GetAttribute("gene_name", value)) {
        return value;
    } else if (GetAttribute("Name", value)) {
        return value;
    } else {
        return gene_id;
    }
}

unsigned GTFFeature::Length() const {
    return (end-start)+1;
}

void GTFFeature::Print() const {
    printf("%s\t%s\t%s\t%u\t%u\t%s\t%c\t%c\t%s", chr.c_str(), source.c_str(), feature.c_str(), start, end, score.c_str(), strand, frame, gene_id.c_str());
    for (id_set::iterator it = transcript_ids.begin(); it != transcript_ids.end(); ++it) {
        printf("\t%s",it->c_str());
    }
    printf("\t[%u]\n", read_count);
}

GTFFeature::GTFFeature(const GTFFeature& rhs) 
    :   key(rhs.key), chr(rhs.chr), source(rhs.source), feature(rhs.feature), start(rhs.start), end(rhs.end), score(rhs.score),
        strand(rhs.strand), frame(rhs.frame), type(rhs.type), gene_id(rhs.gene_id), transcript_id(rhs.transcript_id), 
        gene_name(rhs.gene_name), transcript_name(rhs.transcript_name), transcript_ids(rhs.transcript_ids), 
        attributes(rhs.attributes), read_count(rhs.read_count)
{
    InitializeExclusiveLock(&mutex);
}
  
GTFFeature::~GTFFeature() {
    DestroyExclusiveLock(&mutex);
}

GTFFeature& GTFFeature::operator=(const GTFFeature& rhs) {

    if (this != &rhs) {
        key = rhs.key;
        chr = rhs.chr;
        source = rhs.source;
        feature = rhs.feature;
        start = rhs.start;
        end = rhs.end;
        score = rhs.score;
        strand = rhs.strand;
        frame = rhs.frame;
        type = rhs.type;
        gene_id = rhs.gene_id;
        transcript_id = rhs.transcript_id;
        gene_name = rhs.gene_name;
        transcript_name = rhs.transcript_name;
        transcript_ids = rhs.transcript_ids;
        attributes = rhs.attributes;
        read_count = rhs.read_count;
        InitializeExclusiveLock(&mutex);
    }
    return *this;
}

bool GTFFeature::operator<(const GTFFeature& rhs) const {
    return start < rhs.start;
}

void GTFFeature::IncrementReadCount() {
 
     //Lock the mutex, because multiple threads might try to do this simultaneously
     AcquireExclusiveLock(&mutex);

     read_count++;
 
     //Unlock it
     ReleaseExclusiveLock(&mutex);
 
}

GTFGene::GTFGene(string _chr, string _gene_id, unsigned _start, unsigned _end, string gene_name_) 
    : chr(_chr), gene_id(_gene_id), start(_start), end(_end), gene_name(gene_name_), read_count(0)
{
    InitializeExclusiveLock(&mutex);
}

GTFGene::GTFGene(const GTFGene& rhs) 
    : chr(rhs.chr), gene_id(rhs.gene_id), start(rhs.start), end(rhs.end), gene_name(rhs.gene_name), features(rhs.features), read_count(rhs.read_count)
{
    InitializeExclusiveLock(&mutex);
}
  
GTFGene::~GTFGene() {
    DestroyExclusiveLock(&mutex);
}

GTFGene& GTFGene::operator=(const GTFGene& rhs) {

    if (this != &rhs) {
        chr = rhs.chr;
        gene_id = rhs.gene_id;
        start = rhs.start;
        end = rhs.end;
        gene_name = rhs.gene_name;
        features = rhs.features;
        read_count = rhs.read_count;
        InitializeExclusiveLock(&mutex);
    }
    return *this;
}

bool GTFGene::operator<(const GTFGene &rhs) const {
    printf("< operator not implemented\n");
    exit(0);
}

void GTFGene::UpdateBoundaries(unsigned new_start, unsigned new_end) {
    start = min(new_start, start);
    end = max(new_end, end);
}

void GTFGene::Process(feature_map &all_features, transcript_map &all_transcripts) {

    //Copy all exons into separate vector
    for (id_set::iterator it = transcript_ids.begin(); it != transcript_ids.end(); ++it) {
        
        //Process each transcript
        transcript_map::iterator transcript = all_transcripts.find(*it);
        if (transcript == all_transcripts.end()) {
            printf("Transcript %s not found in complete set", it->c_str());
            exit(1);
        }

        //Pass in all features and my features
        transcript->second.Process(all_features, features);

    }

    /*
    //Print the transcripts
    for (id_set::iterator it = transcript_ids.begin(); it != transcript_ids.end(); ++it) {
        //Process each transcript
        transcript_map::iterator transcript = all_transcripts.find(*it);
        if (transcript == all_transcripts.end()) {
            printf("Transcript %s not found in complete set", it->c_str());
            exit(1);
        }

        //Pass in all features and my features
        transcript->second.Print();
    }
    */
}   

bool GTFGene::CheckBoundary(string query_chr, unsigned query_pos, unsigned buffer) const {

    //Compare chromosomes
    if (chr.compare(query_chr) != 0) {
        return false;
    }
  
    //Compare position within buffered gene coordinates
    if ((query_pos >= max(start-buffer+1, (unsigned)1)) && (query_pos <= end+buffer)) {
        return true;
    }
    return false;
}

void GTFGene::IncrementReadCount() {
    
    //Lock the mutex, because multiple threads might try to do this simultaneously
    AcquireExclusiveLock(&mutex);
        
    read_count++;
        
    //Unlock it
    ReleaseExclusiveLock(&mutex);

}

void GTFGene::WriteJunctionCountID(ofstream &outfile) const {

    for (feature_pointer_map::const_iterator it = features.begin(); it != features.end(); ++it) {
        if (it->second->Type() == INTRON) {
            float gene_expression = ((float)read_count / (float) 1000.0) + 1;
            outfile << gene_id+":"+it->second->Chr()+':'+ToString(it->second->Start())+"-"+ToString(it->second->End()) << '\t' << round((float)it->second->ReadCount() / gene_expression) << endl;
        }
    }
}

void GTFGene::WriteReadCountID(ofstream &outfile) const {
    outfile << gene_id << '\t' << read_count << endl;
}

void GTFGene::WriteReadCountName(ofstream &outfile) const {
    outfile << gene_name << '\t' << read_count << endl;
}

void GTFGene::Print() const {
    printf("%s\t%u\t%u\t%s\n", chr.c_str(), start, end, gene_id.c_str());
}

GTFTranscript::GTFTranscript(string _chr, string _gene_id, string _transcript_id, string _gene_name, string _transcript_name, unsigned _start, unsigned _end) 
    : chr(_chr), gene_id(_gene_id), transcript_id(_transcript_id), gene_name(_gene_name), transcript_name(_transcript_name), start(_start), end(_end), read_count(0)
{
    InitializeExclusiveLock(&mutex);
}

GTFTranscript::GTFTranscript(const GTFTranscript& rhs) 
    : chr(rhs.chr), gene_id(rhs.gene_id), transcript_id(rhs.transcript_id), gene_name(rhs.gene_name), transcript_name(rhs.transcript_name), features(rhs.features), start(rhs.start), end(rhs.end), read_count(rhs.read_count)
{
    InitializeExclusiveLock(&mutex);
}
  
GTFTranscript::~GTFTranscript() 
{
    DestroyExclusiveLock(&mutex);
}

GTFTranscript& GTFTranscript::operator=(const GTFTranscript& rhs) {

    if (this != &rhs) {
        chr = rhs.chr;
        gene_id = rhs.gene_id;
        transcript_id = rhs.transcript_id;
        gene_name = rhs.gene_name;
        transcript_name = rhs.transcript_name;
        features = rhs.features;   
        start = rhs.start;
        end = rhs.end;
        read_count = rhs.read_count;
        InitializeExclusiveLock(&mutex);
    }
    return *this;
}

void GTFTranscript::Process(feature_map &all_features, feature_pointer_map &gene_features) {

    //Copy all exons into separate vector
    sort(features.begin(), features.end(), FeatureListSort);

    feature_list::iterator prev = features.end();
    for (feature_list::iterator current = features.begin(); current != features.end(); ++current) {
        
        //If current is an exon, add it in
        if ((*current)->type == EXON) {

            //If prev is set, also add an intron
            if (prev != features.end()) {

                //Create the new intron with the appropriate boundaries
                GTFFeature intron(**current);
                intron.feature = "intron";
                intron.start = (*prev)->end+1;
                intron.end = (*current)->start-1;
                intron.key = intron.chr + ToString(intron.start) + ToString(intron.end); 
                intron.type = INTRON;

                //Try to add the intron
                feature_map::iterator fpos = all_features.insert(all_features.begin(), feature_map::value_type(intron.key, intron));       
 
                //Insert this transcript_id into this intron
                fpos->second.transcript_ids.insert(intron.transcript_id);

                feature_pointer_map::iterator gpos = gene_features.insert(gene_features.begin(), feature_pointer_map::value_type(intron.key, &fpos->second));

                //Insert this transcript_id into this intron
                gpos->second->transcript_ids.insert(intron.transcript_id);
               
                //Add point of this to the list of exons in this transcript
                exons.push_back(&fpos->second);               
             }
            
           //Finally, add the current exon to the list of exons
           exons.push_back(*current);

           //If this is an exon, then set prev to be current
           prev = current;
        }
    }

    //Sort the exons by start and end
    //sort(exons.begin(), exons.end(), FeatureListSort);
}

/*
void GTFTranscript::Process(feature_map &all_features, feature_list &gene_features) {

    //Copy all exons into separate vector
    for (feature_list::iterator it = features.begin(); it != features.end(); ++it) {
        
        if ((*it)->feature.compare("exon") == 0) {
            exons.push_back(*it);
        }
    }
    
    //Sort the exons by start and end
    sort(exons.begin(), exons.end(), FeatureListSort);

}
*/

void GTFTranscript::UpdateBoundaries(unsigned new_start, unsigned new_end) {
    start = min(new_start, start);
    end = max(new_end, end);
}

void GTFTranscript::Print() const {
    printf("Transcript: %s\n", transcript_id.c_str());
    for (feature_list::const_iterator it = exons.begin(); it != exons.end(); ++it) {
        (*it)->Print();
    }
}   

unsigned GTFTranscript::SplicedLength() const {

    unsigned length = 0;
    for (feature_list::const_iterator it = exons.begin(); it != exons.end(); ++it) {

        //If transcript_pos is less than or equal to this
        length += (*it)->Length();

    }
    //Don't ever return a length of zero
    return max(length, (unsigned)1);
}

void GTFTranscript::IncrementReadCount(unsigned numPotentialTranscripts = 1) {
    
    //Lock the mutex, because multiple threads might try to do this simultaneously
    AcquireExclusiveLock(&mutex);
        
    read_count+= 1.f/(float)numPotentialTranscripts;
        
    //Unlock it
    ReleaseExclusiveLock(&mutex);

}

unsigned GTFTranscript::GenomicPosition(unsigned transcript_pos, unsigned span) const {
        
    //This assumes 1-offset transcript pos, and returns 1-offset genomic position
    //Converts transcript coordinates into genomic coordinates
    for (feature_list::const_iterator it = exons.begin(); it != exons.end(); ++it) {
 
        //If this is an exon
        if ((*it)->type != EXON) {
            continue; 
        }
               
        //If transcript_pos is less than or equal to this
        if (transcript_pos > (*it)->Length()) {
            transcript_pos -= (*it)->Length();
        } else {
        
            //Get current position for start of alignment
            unsigned genome_pos = (*it)->start + transcript_pos - 1;
                
            //Check to see if read exceeds length of transcript.  If so, return 0
            //This is possible due to the way SNAP uses consecutive chromosomes.
            if (genome_pos+span > end) {
                //printf("genome_pos: %u span: %u end: %u\n", genome_pos, span, end);
                //printf("Warning, transcript_pos exceeds transcript length\n");
                return 0;
                
            }
            return genome_pos;
        }
    }
    //printf("Warning, transcript_pos exceeds transcript length\n");
    return 0;
}

void GTFTranscript::Junctions(unsigned transcript_pos, unsigned span, std::vector<junction> &junctions) const {
    
    //Go through each feature of this transcript until we find the start
    unsigned current_pos = 0;
    unsigned end_pos = transcript_pos + span;
    for (feature_list::const_iterator current = exons.begin(); current != exons.end(); ++current) {

        //Get the end of this exon
        if ((*current)->type == EXON) {
            current_pos += (*current)->Length();
        }

        //If we have found the beginning of the transcript
        if (transcript_pos <= current_pos) {

            //If we have crossed an intron, add it to the list of junctions
            if ((*current)->type == INTRON) {
                junctions.push_back(junction(current_pos+1, *current));

            } else if ((*current)->type == EXON) {

                //If transcript_pos is less than or equal to current, the start
                //lies within this feature
                if (current_pos >= end_pos) {
                    return;
                }
            }
        } 
    }
}

/*
void GTFTranscript::Junctions(unsigned transcript_pos, unsigned span, std::vector<junction> &junctions) const {
    
    //printf("Transcript: %s %u %u\n", transcript_id.c_str(), transcript_pos, span);
    
    //Go through each feature of this transcript until we find the start
    unsigned current_pos = 0;
    feature_list::const_iterator current = exons.begin();
    for (feature_list::const_iterator next = ++(exons.begin()); next != exons.end(); ++next) {

        //Get the end of this exon
        current_pos += (*current)->Length();
        
        //printf("Transcript Pos: %u Current Pos: %u \n", transcript_pos,current_pos);
        
        //If transcript_pos is less than or equal to current, the start
        //lies within this feature
        if (transcript_pos <= current_pos) {

            //printf("Inside: %u - %u + 1 >= %u\n", current_pos, transcript_pos, span);

            //Check to see if the span exceeds the current feature. If so, 
            //simply return the current list of junctions
            if (current_pos-transcript_pos+1 >= span) {
                return;
                
            //Otherwise, add the junction to the list of junctions
            } else {
            
                junctions.push_back(junction(current_pos+1, (*next)->start - (*current)->end-1));
                span -= (current_pos-transcript_pos+1);
                transcript_pos += (current_pos-transcript_pos+1);  
            
            }
        }
        current = next;         
    }

}
*/

void GTFTranscript::WriteFASTA(const Genome *genome, std::ofstream &outfile) const {

    //Get the offset for this chromosome
    unsigned offset;
    bool isValid = genome->getOffsetOfContig(chr.c_str(), &offset);

    if (!isValid) {
        printf("Warning: chromosome %s from the annotation is not found in the genome file\n", chr.c_str());
        return;
    }
    
    string sequence;
    for (feature_list::const_iterator it = exons.begin(); it != exons.end(); ++it) {

        if ((*it)->type == EXON) {
    
            //Get the sequence of this feature from the genome
            const char* seq = genome->getSubstring((*it)->start+offset-1,  (*it)->Length());
        
            if (seq == NULL) {
                printf("Warning: transcript %s exceeds boundaries of chromosome %s. Truncating\n", transcript_id.c_str(), chr.c_str());
                //const char* seq = genome->getSubstring((*it)->start+offset-1, (*it)->Length()-amountExceeded, amountExceeded);
                exit(1);
            }
        
            //Temp copy into string variable
            string temp(seq, (*it)->Length());
            sequence += temp;
        }
    }
    outfile << ">" + transcript_id << endl;
    outfile << sequence << endl;

}

unsigned GTFTranscript::NormalizedCount() const {

    return round((float)read_count / ((float)SplicedLength() / 1000.0));

}

void GTFTranscript::WriteReadCountID(ofstream &outfile) const {
    outfile << transcript_id << '\t' << round(read_count) << endl;
}

void GTFTranscript::WriteReadCountName(ofstream &outfile) const {
    outfile << transcript_name << '\t' << round(read_count) << endl;
}

//Constructor
GTFReader::GTFReader(const char* output) {

    prefix = "default";
    if (NULL != output) {
        prefix = output;
        size_t pos;
        if ((pos = prefix.rfind('.')) != std::string::npos) {
            prefix = prefix.substr(0, pos);
        }
    }

}

//Destructor
GTFReader::~GTFReader() {}

void GTFReader::Load(string _filename) {

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

    //Sort each gene, which in turn proceses each transcript
    for (gene_map::iterator it = genes.begin(); it != genes.end(); ++it) {
        it->second.Process(features, transcripts);
    }
    
    //Add all features to interval tree
    for (feature_map::iterator it = features.begin(); it != features.end(); ++it) {
        feature_intervals.push_back(Interval<GTFFeature*>(it->second.start, it->second.end, &it->second));
    }
    feature_tree = IntervalTree<GTFFeature*>(feature_intervals);

    //Add all genes to interval tree
    for (gene_map::iterator it = genes.begin(); it != genes.end(); ++it) {
        gene_intervals.push_back(Interval<GTFGene*>(it->second.start, it->second.end, &it->second));
    }
    gene_tree = IntervalTree<GTFGene*>(gene_intervals);
    
    //Add all transcripts to interval tree
    for (transcript_map::iterator it = transcripts.begin(); it != transcripts.end(); ++it) {
        transcript_intervals.push_back(Interval<GTFTranscript*>(it->second.start, it->second.end, &it->second));
    }
    transcript_tree = IntervalTree<GTFTranscript*>(transcript_intervals);    
    
    _int64 loadTime = timeInMillis() - loadStart;
    printf("%llds. %u features, %u transcripts, %u genes\n", loadTime / 1000, features.size(), transcripts.size(), genes.size());
    
}

int GTFReader::Parse(string line) {

    //Check to see if this is a comment
    if (line[0] == '#') {
        return 1;
    }

    // Create a new GTFFeature from this line on the heap, because we want it to be constant
    GTFFeature feature(line);
    
    if (feature.feature.compare("exon") != 0) {
        return 1;
    }
    
    string value;
    if (!feature.GetAttribute("gene_id", value) && !feature.GetAttribute("Parent", value)) {
        printf("Warning: annotation file missing 'gene_id' (GTF) or 'Parent' (GFF3) for exon entry\n");
    }
    
    //Try to find this feature in the feature map
    //feature_map::iterator fpos = features.insert(feature.key);
    feature_map::iterator fpos = features.insert(features.begin(), feature_map::value_type(feature.key, feature));
    
    //Insert this transcript_id into this feature
    fpos->second.transcript_ids.insert(feature.transcript_id);

    //Try to find this transcript in the transcript_map
    transcript_map::iterator pos = transcripts.find(feature.transcript_id);
        
    //If this sequence is not found, create a new vector to store this sequence (and others like it)
    if ((pos == transcripts.end())) {
        GTFTranscript transcript(feature.chr, feature.gene_id, feature.transcript_id, feature.GeneName(), feature.TranscriptName(), feature.start, feature.end);
        transcript.features.push_back(&fpos->second);
        transcripts.insert(transcript_map::value_type(feature.transcript_id, transcript));
        
    //Otherwise, add this feature to the transcript
    } else {
        pos->second.features.push_back(&fpos->second);
        pos->second.UpdateBoundaries(feature.start, feature.end);
    }
    
    //Try to find this gene in the gene_map
    gene_map::iterator gpos = genes.find(feature.gene_id);
    
    //If this sequence is not found, create a new vector to store this sequence (and others like it)
    if ((gpos == genes.end())) {
        GTFGene gene(feature.chr, feature.gene_id, feature.start, feature.end, feature.GeneName());
        //gene.features.push_back(&fpos->second);
        gene.transcript_ids.insert(feature.transcript_id);
        genes.insert(gene_map::value_type(feature.gene_id, gene));
        
    //Otherwise, add this feature to the gene
    } else {
        //gpos->second.features.push_back(&fpos->second);
        gpos->second.transcript_ids.insert(feature.transcript_id);
        gpos->second.UpdateBoundaries(feature.start, feature.end);
    }    
       
    return 0;
    
}

const GTFTranscript& GTFReader::GetTranscript(string transcript_id) const {
 
    transcript_map::const_iterator pos = transcripts.find(transcript_id);
    if (pos == transcripts.end()) {
        //raise exception
        printf("No transcript %s\n", transcript_id.c_str());
        exit(1);
    }
    return pos->second;

}

const GTFGene& GTFReader::GetGene(string gene_id) const {
 
    gene_map::const_iterator pos = genes.find(gene_id);
    if (pos == genes.end()) {
        //raise exception
        printf("No gene %s\n", gene_id.c_str());
        exit(1);
    }
    return pos->second;

}

void GTFReader::IncrementReadCount(string transcript_id0, unsigned transcript_start0, unsigned start0, unsigned length0) {

    transcript_map::iterator pos = transcripts.find(transcript_id0);
    if (pos == transcripts.end()) {
        //raise exception
        printf("No transcript %s\n", transcript_id0.c_str());
        exit(1);
    }
    string gene_id = pos->second.GeneID();

    //Increment the gene count for one of the transcripts
    gene_map::iterator gpos = genes.find(gene_id);
    if (gpos == genes.end()) {
        //raise exception
        printf("No gene %s\n", gene_id.c_str());
        exit(1);
    }
    gpos->second.IncrementReadCount();

}

void GTFReader::IncrementReadCount(string transcript_id0, unsigned transcript_start0, unsigned start0, unsigned length0, string transcript_id1, unsigned transcript_start1, unsigned start1, unsigned length1) {

    //If a read is aligned to transcriptome, it can be spliced
    //If it is not aligned to transcriptome, it cannot be spliced
    std::set<string> transcript_ids0, transcript_ids1;
    std::set<string>::iterator pos;
    
    //If the first read is aligned to transcriptome
    if (transcript_id0.size() != 0) {
    
        //Get the transcript associated with this id
        const GTFTranscript &transcript0 = GetTranscript(transcript_id0);
        
        //Get the junctions for this transcript
        std::vector<junction> junctions;
        transcript0.Junctions(transcript_start0, length0, junctions);
        
        for (std::vector<junction>::iterator it = junctions.begin(); it != junctions.end(); ++it) {
 
            //Increment the splice count for this junction
            it->second->IncrementReadCount();
       
            unsigned length = it->first - transcript_start0;
            
            //printf("Querying with [%u %u]\n", start0, start0+length-1);
            
            //For each junction, query the feature tree
            std::vector<GTFFeature*> results;
            IntervalFeatures(transcript0.chr, start0, start0+length-1, results);
                                                
            if (transcript_ids0.size() == 0) {
            
                for (std::vector<GTFFeature*>::iterator it2 = results.begin(); it2 != results.end(); ++it2) {
                    transcript_ids0.insert((*it2)->transcript_id);
                }               
                 
            } else {
                                                
                //For each result, add the transcript_id to the set only if it already exists
                std::set<string> temp_ids;
                for (std::vector<GTFFeature*>::iterator it2 = results.begin(); it2 != results.end(); ++it2) {
                    if ((pos = transcript_ids0.find((*it2)->transcript_id)) != transcript_ids0.end()) {
                        temp_ids.insert((*it2)->transcript_id);
                    }
                }
                transcript_ids0 = temp_ids;
                
            }
                        
            transcript_start0 += length;
            start0 += length + it->second->Length();
            length0 -= length; 
             
        } 
        
        //Do the last junction (or if there were no junctions, do the entire read)
        //printf("Querying with [%u %u]\n", start0, start0+length0-1);
        std::vector<GTFFeature*> results;
        IntervalFeatures(transcript0.chr, start0, start0+length0-1, results);

        if (transcript_ids0.size() == 0) {
            for (std::vector<GTFFeature*>::iterator it2 = results.begin(); it2 != results.end(); ++it2) {
                transcript_ids0.insert((*it2)->transcript_id);
            }               
             
        } else {
                                            
            //For each result, add the transcript_id to the set only if it already exists
            std::set<string> temp_ids;
            for (std::vector<GTFFeature*>::iterator it2 = results.begin(); it2 != results.end(); ++it2) {
                if ((pos = transcript_ids0.find((*it2)->transcript_id)) != transcript_ids0.end()) {
                    temp_ids.insert((*it2)->transcript_id);
                }
            }
            transcript_ids0 = temp_ids;
            
        }

    } else {
        //printf("You have not implemented if the read is aligned to the genome!\n");
        return;
    }
        
    //If the second read is aligned to transcriptome
    if (transcript_id1.size() != 0) {
    
        //Get the transcript associated with this id
        const GTFTranscript &transcript1 = GetTranscript(transcript_id1);
        
        //Get the junctions for this transcript
        std::vector<junction> junctions;
        transcript1.Junctions(transcript_start1, length1, junctions);
        
        for (std::vector<junction>::iterator it = junctions.begin(); it != junctions.end(); ++it) {

            //Increment the splice count for this junction           
            it->second->IncrementReadCount();

            unsigned length = it->first - transcript_start1;
            
            //printf("Querying with [%u %u]\n", start1, start1+length-1);
            
            //For each junction, query the feature tree
            std::vector<GTFFeature*> results;
            IntervalFeatures(transcript1.chr, start1, start1+length-1, results);
            
            if (transcript_ids1.size() == 0) {
                for (std::vector<GTFFeature*>::iterator it2 = results.begin(); it2 != results.end(); ++it2) {
                    transcript_ids1.insert((*it2)->transcript_id);
                }               
                 
            } else {
                                                
                //For each result, add the transcript_id to the set only if it already exists
                std::set<string> temp_ids;
                for (std::vector<GTFFeature*>::iterator it2 = results.begin(); it2 != results.end(); ++it2) {
                    if ((pos = transcript_ids1.find((*it2)->transcript_id)) != transcript_ids1.end()) {
                        temp_ids.insert((*it2)->transcript_id);
                    }
                }
                transcript_ids1 = temp_ids;
                
            }
            
            transcript_start1 += length;
            start1 += length + it->second->Length();
            length1 -= length; 
            
        } 
        
        //Do the last junction (or if there were no junctions, do the entire read)
        //printf("Querying with [%u %u]\n", start1, start1+length1-1);
        std::vector<GTFFeature*> results;
        IntervalFeatures(transcript1.chr, start1, start1+length1-1, results);
        
        if (transcript_ids1.size() == 0) {
        
            for (std::vector<GTFFeature*>::iterator it2 = results.begin(); it2 != results.end(); ++it2) {
                transcript_ids1.insert((*it2)->transcript_id);
            }               
             
        } else {
                                            
            //For each result, add the transcript_id to the set only if it already exists
            std::set<string> temp_ids;
            for (std::vector<GTFFeature*>::iterator it2 = results.begin(); it2 != results.end(); ++it2) {
                if ((pos = transcript_ids1.find((*it2)->transcript_id)) != transcript_ids1.end()) {
                    temp_ids.insert((*it2)->transcript_id);
                }
            }
            transcript_ids1 = temp_ids;
            
        }
    
    } else {
        //printf("You have not implemented if the read is aligned to the genome!\n");
        return;
    }
    
    //Finally, only increment the counts for those that are defined in both lists
    std::set<string> final_ids;
    for (std::set<string>::iterator it = transcript_ids0.begin(); it != transcript_ids0.end(); ++it) {
        if ((pos = transcript_ids1.find(*it)) != transcript_ids1.end()) {
            final_ids.insert(*it);
        }
    }
    
    if (final_ids.size() == 0) {
        //printf("Warning, no feature discovered for alignment\n");
        //printf("%s %s\n", transcript_id0.c_str(), transcript_id1.c_str());
        return;
    }
    
    string gene_id;
    //Finally, increment the transcript count for these transcript(s)
    for (std::set<string>::iterator it = final_ids.begin(); it != final_ids.end(); ++it) {
        transcript_map::iterator pos = transcripts.find(*it);
        if (pos == transcripts.end()) {
            //raise exception
            printf("No transcript %s\n", it->c_str());
            exit(1);
        }
        gene_id = pos->second.GeneID();
        //We only increment once for a paired-end fragment
        pos->second.IncrementReadCount((unsigned) final_ids.size());
    }
    
    //Increment the gene count for one of the transcripts
    gene_map::iterator gpos = genes.find(gene_id);
    if (gpos == genes.end()) {
        //raise exception
        printf("No gene %s\n", gene_id.c_str());
        exit(1);
    }
    //We only increment once for a paired-end fragment
    gpos->second.IncrementReadCount();    

//     printf("Final\n");
//     for (std::set<string>::iterator it = final_ids.begin(); it != final_ids.end(); ++it) {
//         printf("%s\n", it->c_str());
//     }
    
}


void GTFReader::IntervalFeatures(std::string chr, unsigned start, unsigned stop, std::vector<GTFFeature*> &results) {

    std::vector<Interval<GTFFeature*> > temp_results;
    feature_tree.findOverlapping(start, stop, temp_results);

    //Filter by chromosome
    for (std::vector<Interval<GTFFeature*> >::iterator it = temp_results.begin(); it != temp_results.end(); ++it) {
        if (it->value->chr.compare(chr) == 0) {
            results.push_back(it->value);
        }
    }
}

void GTFReader::IntervalTranscripts(std::string chr, unsigned start, unsigned stop, std::vector<GTFTranscript> &results) {
    
    std::vector<Interval<GTFTranscript*> > temp_results;
    transcript_tree.findOverlapping(start, stop, temp_results);

    //Filter by chromosome
    for (std::vector<Interval<GTFTranscript*> >::iterator it = temp_results.begin(); it != temp_results.end(); ++it) {
        if (it->value->chr.compare(chr) == 0) {
            results.push_back(*(it->value));
        }
    }
}

void GTFReader::IntervalGenes(std::string chr, unsigned start, unsigned stop, std::vector<GTFGene> &results) {
    
    std::vector<Interval<GTFGene*> > temp_results;
    gene_tree.findOverlapping(start, stop, temp_results);

    //Filter by chromosome
    for (std::vector<Interval<GTFGene*> >::iterator it = temp_results.begin(); it != temp_results.end(); ++it) {
        if (it->value->chr.compare(chr) == 0) {
            results.push_back(*(it->value));
        }
    }
}

void GTFReader::IntrageneUnannotatedPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq0, string seq1) {
    intragene_unannotated_pairs.AddInterval(chr0, start0, end0, chr1, start1, end1, id, seq0, seq1, false);
}

void GTFReader::IntrageneUnannotatedSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq) {
    intragene_unannotated_splices.AddInterval(chr0, start0, end0, chr1, start1, end1, id, seq, "", true);
}

void GTFReader::IntrageneCircularPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq0, string seq1) {
    intragene_circular_pairs.AddInterval(chr0, start0, end0, chr1, start1, end1, id, seq0, seq1, false);
}

void GTFReader::IntrageneCircularSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq) {
    intragene_circular_splices.AddInterval(chr0, start0, end0, chr1, start1, end1, id,  seq, "", true);
}

void GTFReader::IntrachromosomalPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq0, string seq1) {
   intrachromosomal_pairs.AddInterval(chr0, start0, end0, chr1, start1, end1, id, seq0, seq1, false);
}

void GTFReader::IntrachromosomalSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq) {
    intrachromosomal_splices.AddInterval(chr0, start0, end0, chr1, start1, end1, id, seq, "", true);
}

void GTFReader::InterchromosomalPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq0, string seq1) {
    interchromosomal_pairs.AddInterval(chr0, start0, end0, chr1, start1, end1, id, seq0, seq1, false);
}

void GTFReader::InterchromosomalSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, string seq) {
    interchromosomal_splices.AddInterval(chr0, start0, end0, chr1, start1, end1, id, seq, "", true);

  /*
  string gene1 = "BCR";
  string gene2 = "ABL1";
 
  std::vector<GTFGene> results0, results1;
   IntervalGenes(chr0, start0, end0, results0);
   IntervalGenes(chr1, start1, end1, results1);
   for (std::vector<GTFGene>::iterator it = results0.begin(); it != results0.end(); ++it) {
     if ((it->gene_name.compare(gene1) == 0) || (it->gene_name.compare(gene2) == 0)) {
       
       for (std::vector<GTFGene>::iterator it2 = results1.begin(); it2 != results1.end(); ++it2) {
         if ((it2->gene_name.compare(gene1) == 0) || (it2->gene_name.compare(gene2) == 0)) {
          found = true; 
          printf("SPLICE Found BCR_ABL\n");
         }
       }
     }
   }
  */
}

void GTFReader::WriteReadCounts() {

    printf("Writing gene and junction read counts... ");
    fflush(stdout);
    _int64 loadStart = timeInMillis();

    ofstream transcript_name_counts, transcript_id_counts;
    ofstream gene_name_counts, gene_id_counts;
    ofstream junction_name_counts, junction_id_counts;

    //Open output files for id counts
    transcript_id_counts.open((prefix+".transcript_id.counts.txt").c_str());
    gene_id_counts.open((prefix+".gene_id.counts.txt").c_str());
    junction_id_counts.open((prefix+".junction_id.counts.txt").c_str());    

    transcript_name_counts.open((prefix+".transcript_name.counts.txt").c_str());
    gene_name_counts.open((prefix+".gene_name.counts.txt").c_str());
    junction_name_counts.open((prefix+".junction_name.counts.txt").c_str());
    
    //Go though each transcript and write
    for (transcript_map::iterator it = transcripts.begin(); it != transcripts.end(); ++it) {
        it->second.WriteReadCountID(transcript_id_counts);
        it->second.WriteReadCountName(transcript_name_counts);
    }
    
    //Go though each transcript and write
    for (gene_map::iterator it = genes.begin(); it != genes.end(); ++it) {
        it->second.WriteReadCountID(gene_id_counts);
        it->second.WriteJunctionCountID(junction_id_counts);
    }
    
    //Create a map to store the gene names
    std::map<string, unsigned> gene_counts;
    std::map<string, unsigned>::iterator pos;
    
    //Go though each gene and write
    for (gene_map::iterator it = genes.begin(); it != genes.end(); ++it) {
        //it->second.WriteReadCountID(gene_id_counts);
        //it->second.WriteReadCountName(gene_name_counts);
        //Get the gene name
        string gene_name = it->second.GeneName();
        
        pos = gene_counts.find(gene_name);
        //If the gene is not found, initialize the read count with the current value
        if (pos == gene_counts.end()) {
            gene_counts.insert(std::map<string, unsigned>::value_type(gene_name, it->second.ReadCount()));
            
        //If the gene is found, increment it by the current gene count
        } else {
            pos->second += it->second.ReadCount();
        }
    }
    
    //Finally, write out the consolidated genes and read counts
    for (pos = gene_counts.begin(); pos != gene_counts.end(); ++pos) {
        gene_name_counts << pos->first << '\t' << pos->second << endl;
    
    }
    
    transcript_id_counts.close();
    gene_id_counts.close();
    junction_id_counts.close();
    transcript_name_counts.close();
    gene_name_counts.close();
    junction_name_counts.close();

    _int64 loadTime = timeInMillis() - loadStart;
    printf("%llds.\n", loadTime / 1000);
    fflush(stdout);
}

void GTFReader::AnalyzeReadIntervals() {

    printf("Calculating gene fusions... ");
    fflush(stdout);
    _int64 loadStart = timeInMillis();

    bool filterPromiscuousGenes = true;
    unsigned pairedBuffer = 100;
    unsigned splicedBuffer = 0;
    
    unsigned minCount = 5;
    unsigned intersectionBuffer = 10;

    //Open output file
    ofstream interchromosomal_file, intrachromosomal_file, unannotated_file, circular_file;
    ofstream logfile, readfile;
    interchromosomal_file.open((prefix+".interchromosomal.fusions.gtf").c_str());
    intrachromosomal_file.open((prefix+".intrachromosomal.fusions.gtf").c_str());
    //unannotated_file.open((prefix+".unannotated_intervals.gtf").c_str());
    //circular_file.open((prefix+".circular_intervals.gtf").c_str());
    logfile.open((prefix+".fusions.txt").c_str());
    readfile.open((prefix+".fusions.reads.fa").c_str());

    interchromosomal_pairs.Consolidate(this, pairedBuffer);
    interchromosomal_splices.Consolidate(this, splicedBuffer);
    interchromosomal_splices.Intersect(interchromosomal_pairs, intersectionBuffer, minCount, this);
    logfile << "Inter-Chromosomal Fusions" << endl;
    interchromosomal_splices.WriteGTF(interchromosomal_file);
    interchromosomal_splices.WriteSplicedMatePairs(logfile, readfile);
    logfile << endl;
    interchromosomal_splices.Clear();
    
    intrachromosomal_pairs.Consolidate(this, pairedBuffer);
    intrachromosomal_splices.Consolidate(this, splicedBuffer);
    intrachromosomal_splices.Intersect(intrachromosomal_pairs, intersectionBuffer, minCount, this);
    logfile << "Intra-Chromosomal Fusions" << endl;
    intrachromosomal_splices.WriteGTF(intrachromosomal_file);
    intrachromosomal_splices.WriteSplicedMatePairs(logfile, readfile);
    logfile << endl;
    intrachromosomal_splices.Clear();
    
    /*
    intragene_unannotated_pairs.Consolidate(this, pairedBuffer);
    intragene_unannotated_splices.Consolidate(this, splicedBuffer);
    intragene_unannotated_pairs.Intersect(intragene_unannotated_splices, intersectionBuffer, minCount, this);
    logfile << "Intra-Gene Unannotated Intervals" << endl;
    intragene_unannotated_pairs.WriteGTF(unannotated_file);
    intragene_unannotated_pairs.WriteSplicedMatePairs(logfile, readfile);
    logfile << endl;
    intragene_unannotated_pairs.Clear();
    
    intragene_circular_pairs.Consolidate(this, pairedBuffer);
    intragene_circular_splices.Consolidate(this, splicedBuffer);
    intragene_circular_pairs.Intersect(intragene_circular_splices, intersectionBuffer, minCount, this); 
    logfile << "Intra-Gene Circular Intervals" << endl;
    intragene_circular_pairs.WriteGTF(circular_file);
    intragene_circular_pairs.WriteSplicedMatePairs(logfile, readfile);
    logfile << endl;
    intragene_circular_pairs.Clear();
    */    

    interchromosomal_file.close();
    intrachromosomal_file.close();
    //unannotated_file.close();
    //circular_file.close();
    logfile.close();
    //readfile.close();
    
    _int64 loadTime = timeInMillis() - loadStart;
    printf("%llds.\n", loadTime / 1000);
    fflush(stdout);
}

void GTFReader::BuildTranscriptome(const Genome *genome) {

    printf("Building transcriptome FASTA file... ");
    fflush(stdout);
    _int64 loadStart = timeInMillis();

    //Save the input filename
    string filename = "transcriptome.fa";
        
    //Open input file
    std::ofstream outfile(filename.c_str(), std::ios::out);
    if (!outfile.is_open()) { 
        printf("Cannot open output file '%s'\n", filename.c_str());
	    exit(1);
    }
    
    //Loop over all transcripts
    for (transcript_map::iterator it = transcripts.begin(); it != transcripts.end(); ++it) { 
        it->second.WriteFASTA(genome, outfile); 
    }
    
    //Close the output file
    outfile.close();

    _int64 loadTime = timeInMillis() - loadStart;
    printf("%llds.\n", loadTime / 1000);

}

/*
int main(int argc, char *argv[]) {

    const char *fasta_filename = NULL;
    const char *gtf_filename = NULL;

    //Process the command-line arguments
    int c;
    while ((c=getopt(argc, argv, "g:f:")) != EOF) {
    
        switch(c) {
            case 'g':
                gtf_filename = optarg;
                break;
            case 'f':
                fasta_filename = optarg;
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
    gtf.Load(gtf_filename);
    
    const Genome *genome = ReadFASTAGenome(fasta_filename);
    //gtf.Test();
    
    struct timeval t2;
    gettimeofday(&t2, NULL);
    
    double elapsedTime;
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;
    fprintf(stderr, "Time elapsed: %.3lf seconds\n", elapsedTime / 1000.0);
    return 0;
}
*/

/*
void GTFReader::Test() {

    
    const GTFTranscript &transcript = GetTranscript("ENST00000489673");

    std::vector<junction> junctions = transcript.Junctions(1, 200);
    
    for (std::vector<junction>::iterator it = junctions.begin(); it != junctions.end(); ++it) {
        printf("[%d %d]\n", it->first, it->second);
    }
    
    std::vector<GTFGene> results;
    IntervalGenes("22", 42290565, 42359490, results);
    printf("Size: %d\n", results.size());
    for (std::vector<GTFGene>::iterator it = results.begin(); it != results.end(); ++it) {
        it->Print();
    }
}
*/

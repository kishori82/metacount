#include "core.h"
#include "types.h"

using namespace std;

#define _MAX 1000000000

CONTIG_ORF *create_contig_orf_map(const string &orf_file) {
    CONTIG_ORF *contig_orf;

    if ((contig_orf = new CONTIG_ORF) == nullptr) {
       error("ERROR: Failed to allocate CONTIG_ORF");
    }

    MatchOutputParser *parser = ParserFactory::createParser(orf_file, std::string("GFF"));

    MATCH  match;
    std::map<string, vector<MATCH> > orf_dictionary;
    vector<MATCH> match_vector;

    // reads in all the orfs and associated to the contigs 
    std::cout << "Sorting through the reads...." ;
    string contigid, orfid;
    for(int i =0; ; i++ )  {
       if( !parser->nextline(match) )  break;
/*
       std::cout << shorten_id(match.query, CONTIGID) 
                 << "\t" <<  match.query 
                 << "\t" << match.subject << "\t" 
                 << "\t" << shorten_id(match.subject, ORFID) << "\t" 
                 << "\t" << shorten_id(match.subject, ORFID) << "\t" 
                 << match.start << "\t"<< match.end << std::endl;
*/

       contigid =  shorten_id(match.query, CONTIGID); 
       orfid =  shorten_id(match.subject, ORFID); 

       if (contig_orf->find(contigid) == contig_orf->end()) {
          contig_orf->insert(std::pair<string, vector<ORFINFO *>*>(contigid, new std::vector<ORFINFO *> ));
       }
//       std::cout << "\t" << orfid << std::endl;

       assert(contig_orf->find(contigid) != contig_orf->end());

       auto it = contig_orf->find(contigid);
       ORFINFO *orfinfo = new ORFINFO(orfid, match.start, match.end);
       it->second->push_back(orfinfo);

/*
       if(i%10000==0) {
           std::cout << "\n\033[F\033[J";
       }
*/
    }
    delete parser;
    return contig_orf;
}


std::map<string, uint32_t> * get_contig_information(vector<string>  read_map_files) {
  SamFile samIn;
  SamFileHeader samFileHeader;
  SamHeaderRecord *headerRec;//= samFileHeader.getNextSQRecord();	
  string contig_name;
  uint32_t contig_len;

  std::map<string, uint32_t> * contig_length = new std::map<string, uint32_t>;

  for (uint32_t i = 0; i < read_map_files.size(); i++) {
     samIn.OpenForRead(read_map_files[i].c_str());
     samIn.ReadHeader(samFileHeader);
     while ((headerRec = samFileHeader.getNextSQRecord()) != NULL) {
       //std::cout << headerRec->getTagValue("SN") << "\t" << headerRec->getTagValue("LN") <<  std::endl;
       contig_name = headerRec->getTagValue("SN");
       contig_len = strToUint16(headerRec->getTagValue("LN"));
       contig_length->insert(std::pair<string, uint32_t>(contig_name, contig_len));       
     }
  }
  return contig_length;
}


ORFINFO *find_orf_match(vector<ORFINFO *> *orf_vector, std::pair<uint32_t, uint32_t> &loc_pair) {
  uint32_t min_loc, max_loc, min_orf, max_orf;

  min_loc = min(loc_pair.first, loc_pair.second);
  max_loc = max(loc_pair.first, loc_pair.second);

  for (auto it = orf_vector->begin(); it != orf_vector->end(); it++) {
    min_orf = min((*it)->start, (*it)->end);
    max_orf = max((*it)->start, (*it)->end);

    if (max_loc < min_orf) {
       return nullptr;
    }

    if (max_loc >= min_orf && min_loc <= max_orf) {
       return *it;
    } 
  }
  return nullptr;
}

std::map<string, uint32_t> * read_bam_files(CONTIG_ORF *contig_orf, 
                                            vector<string> read_map_files)
{
  SamFile samIn;
  SamFileHeader samFileHeader;
  SamRecord samRecord;

  string contigid;
  std::pair<uint32_t, uint32_t> loc_pair;

  ORFINFO *destorf;
  vector<ORFINFO *> *orf_vector;

  MATCH  match;
  std::map<string, uint32_t> * contig_read_counts = new std::map<string, uint32_t>;
  for (size_t i = 0; i < read_map_files.size(); i++) {
     std::cout << "processing : " << read_map_files[i] << std::endl;
     samIn.OpenForRead(read_map_files[i].c_str());
     samIn.ReadHeader(samFileHeader);

     RUN_STATS *sam_file_stats  = new RUN_STATS;
     while (samIn.ReadRecord(samFileHeader, samRecord)) {
       contigid =  shorten_id(samRecord.getReferenceName(), CONTIGID); 

       sam_file_stats->num_total_alns++;
       
       if (contig_read_counts->find(contigid) ==  contig_read_counts->end()) {
         contig_read_counts->insert(std::pair<string, uint32_t>(contigid, 0));
       }

       (contig_read_counts->find(contigid)->second)++;

       // get the information for the the alignment 
       // these four immediate attributes  are not used,  can be remove
       // useful for printing out information 
       match.query = samRecord.getReadName();
       match.subject = samRecord.getReferenceName();
       match.start = samRecord.get1BasedPosition();
       match.end = samRecord.get1BasedAlignmentEnd();

       // these are required
       match.alnlength = samRecord.getAlignmentLength();
       getSamFlagInfo(samRecord.getFlag(), match);

       if (match.readtype == SE && match.chimeric == NON_SUPPLEMENTARY) { // & (match.aligntype == PRIMARY)  ) {
          sam_file_stats->num_reads_1++;
       } else {
           //  std::cout << match.print() << std::endl;
       }

       if (match.readtype == PE && match.readno == FIRST && match.chimeric == NON_SUPPLEMENTARY) { // & (match.aligntype == PRIMARY)  ) {
          sam_file_stats->num_reads_1++;
       } else {
          if (match.readtype == PE && match.readno == FIRST) {
       //      std::cout << match.print() << std::endl;
          }
       }
     
       if (match.readno == SECOND && match.chimeric == NON_SUPPLEMENTARY) { // & (match.aligntype == PRIMARY)  ) {
          sam_file_stats->num_reads_2++;
       } else {
          //std::cout << match.print() << std::endl;
       }

       if (match.alnvalid == VALID) {
          sam_file_stats->num_mapped_alns++;
       } else {
          sam_file_stats->num_unmapped_alns++;
       }
       

       // if the ORF is in the contig-to-orf map then add the count
       if (contig_orf->find(contigid) !=  contig_orf->end()) {
         loc_pair.first  =  samRecord.get1BasedPosition();
         loc_pair.second =  samRecord.get1BasedAlignmentEnd();

         // if the ORF is in the contig-to-orf map then add the count
          orf_vector = contig_orf->find(contigid)->second;
          destorf = find_orf_match(orf_vector, loc_pair);

         // if found a valid overlapping contig then add the count
         if (destorf != nullptr) {
            destorf->count = destorf->count + match.weight;
         }
       }
      
#ifdef PRINT_VERBOSE
       std::cout << samRecord.getReferenceName() 
                 << "\t" << samRecord.get1BasedPosition()
                 << "\t" << samRecord.get1BasedAlignmentEnd()
                 << "\t" << samRecord.getAlignmentLength() 
                 << "\t" << samRecord.get1BasedPosition() + samRecord.getAlignmentLength() - 1
                 << std::endl;
#endif
     }

     sam_file_stats->print_stats(&std::cout);
  }

  return contig_read_counts;
}

std::vector<TRIPLET>::iterator binary_search(std::vector<TRIPLET> &A, int seekValue)
{
  // continually narrow search until just one element remains
  unsigned int imin, imax;
  imin = 0; imax = A.size();

  while (imin < imax)
    {
      unsigned int imid = (imin+imax)/2;
 
      // code must guarantee the interval is reduced at each iteration
      assert(imid < imax);

      // note: 0 <= imin < imax implies imid will always be less than imax
 
      // reduce the search
      if (A[imid].start < static_cast<unsigned int>(seekValue) )
        imin = imid + 1;
      else
        imax = imid;
    }

    std::vector<TRIPLET>::iterator it = A.begin() + imin;
 //   std::cout <<  it->start << "   " << seekValue << std::endl;
    return it ;
}

void writeOut_ORFwise_RPKM_values(const string orf_rpkm_file,  map<string, float> &orfnames) {
    ofstream output_file;

    output_file.open(orf_rpkm_file.c_str());
    output_file << "# ORF_ID\tCOUNT" << endl;
    for(map<string, float>::iterator it= orfnames.begin(); it != orfnames.end(); it++)  {
         output_file << it->first << "\t" << it->second << endl;
    }
    output_file.close();

}

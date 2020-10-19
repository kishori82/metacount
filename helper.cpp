#include "helper.h"
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


std::map<string, uint16_t> * get_contig_information(vector<string>  read_map_files) {
  SamFile samIn;
  SamFileHeader samFileHeader;
  SamHeaderRecord *headerRec;//= samFileHeader.getNextSQRecord();	
  string contig_name;
  uint16_t contig_len;

  std::map<string, uint16_t> * contig_length = new std::map<string, uint16_t>;

  for (uint32_t i = 0; i < read_map_files.size(); i++) {
     samIn.OpenForRead(read_map_files[i].c_str());
     samIn.ReadHeader(samFileHeader);
     while ((headerRec = samFileHeader.getNextSQRecord()) != NULL) {
       //std::cout << headerRec->getTagValue("SN") << "\t" << headerRec->getTagValue("LN") <<  std::endl;
       contig_name = headerRec->getTagValue("SN");
       contig_len = strToUint16(headerRec->getTagValue("LN"));
       contig_length->insert(std::pair<string, uint16_t>(contig_name, contig_len));       
     }
     std::cout << " get num SQs " << samFileHeader.getNumSQs() << std::endl;
  }
  return contig_length;
}


ORFINFO *find_orf_match(vector<ORFINFO *> *orf_vector, std::pair<uint16_t, uint16_t> &loc_pair) {
  uint16_t min_loc, max_loc, min_orf, max_orf;

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
  std::pair<uint16_t, uint16_t> loc_pair;

  ORFINFO *destorf;
  vector<ORFINFO *> *orf_vector;

  std::map<string, uint32_t> * contig_read_counts = new std::map<string, uint32_t>;

  for (size_t i = 0; i < read_map_files.size(); i++) {
     samIn.OpenForRead(read_map_files[i].c_str());
     samIn.ReadHeader(samFileHeader);
     while (samIn.ReadRecord(samFileHeader, samRecord)) {
       contigid =  shorten_id(samRecord.getReferenceName(), CONTIGID); 

       if (contig_read_counts->find(contigid) ==  contig_read_counts->end()) {
         contig_read_counts->insert(std::pair<string, uint32_t>(contigid, 0));
       }
       (contig_read_counts->find(contigid)->second)++;

       
       if (contig_orf->find(contigid) !=  contig_orf->end()) {
         loc_pair.first  =  samRecord.get1BasedPosition();
         loc_pair.second =  samRecord.get1BasedAlignmentEnd();

         orf_vector = contig_orf->find(contigid)->second;
         destorf = find_orf_match(orf_vector, loc_pair);
         if (destorf != nullptr) {
           destorf->count++;
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
  }

  return contig_read_counts;
}


unsigned int  getMaxReadSize(
       std::map<string, vector<MATCH> > &orf_dictionary,
       std::map<string, CONTIG> &contigs_dictionary
    ) 
    {
    unsigned int size = 0;
    std::map<string, vector<MATCH> >::iterator itcont;

    for(itcont= orf_dictionary.begin(); itcont != orf_dictionary.end(); itcont++)  {
       for(std::vector<TRIPLET>::iterator it = contigs_dictionary[itcont->first].M.begin(); it != contigs_dictionary[itcont->first].M.end(); it++) {
          if( size < it->end  - it->start ) size = it->end  - it->start ;
       }
    }

    return size;
} 


std::vector<TRIPLET>::iterator  binary_search(std::vector<TRIPLET> &A, int seekValue)
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

void substring_coverage(std::map<string, CONTIG> &contigs_dictionary, const std::string &contig,\
                          unsigned long start, unsigned long end, COVERAGE &coverage, unsigned int maxReadLength) {


      if( contigs_dictionary.find(contig) == contigs_dictionary.end()   || contigs_dictionary[contig].L==0 ) {
          coverage.coverage =0 ;
          coverage.numreads =0 ;
          coverage.substring_length = end - start ;
          coverage.uncovered_length = 0;
      }

      float numreads =0; float uncovered_length = 0; float _coverage = 0;
      unsigned long p_end = start;

      int _seekValue =  maxReadLength == 0 ? 0 :  (start < maxReadLength ? 0 : start-maxReadLength );


     std::vector<TRIPLET>::iterator it= contigs_dictionary[contig].M.begin(); 

     if( _seekValue >0 )
        it =  binary_search(contigs_dictionary[contig].M, _seekValue);

      
     //iterate through every read in that contig
     for( ; it != contigs_dictionary[contig].M.end(); it++) {
         uncovered_length  +=  ( p_end > it->start  || it->start > end) ? 0 : it->start - p_end;
        //make sure the read start and end are not going beyoind the contig
         if( it->end > p_end ) p_end = it->end;  
          
         if( (start <= it->start && it->start <= end) ||  (start <= it->end && it->end <= end)  ) {
            numreads += 1;
            // KK originalnumreads += 1/static_cast<float>(it->multi);
      //      std::cout << "multireads " << std::endl;
         }

         // the subsequent reads are going off the end of the orf
         if( it->start > end ) break;

     }
     uncovered_length += (p_end > end ) ? 0 : (end - p_end);

     unsigned long substring_length = end - start;
     if( substring_length > 0 ) 
        _coverage = ((float)(substring_length - uncovered_length )/(float)substring_length)*100;

     coverage.numreads = numreads;
     coverage.coverage = _coverage;
     coverage.substring_length  = substring_length;
     coverage.uncovered_length =  uncovered_length;
}


unsigned long ORFWise_coverage(map<string, CONTIG> &contigs_dictionary, const string &orf_file,\
                               map<string, float> &orfnames, unsigned long genome_length,\
                               unsigned long &orf_length,  unsigned long num_mappable_reads,\
                               int count_type, float genome_equivalent,  bool show_status) {

    MATCH  match;
    COVERAGE coverage;
    MatchOutputParser *parser = ParserFactory::createParser(orf_file, std::string("GFF"));

    std::map<string, vector<MATCH> > orf_dictionary;

    unsigned long _num_orfs=0;

    vector<MATCH> match_vector;

    // reads in all the orfs and associated to the contigs 
    std::cout << "Sorting through the reads...." ;
    for(int i =0; ; i++ )  {
       if( !parser->nextline(match) )  break;
       if(show_status && i%10000==0) {
           std::cout << "\n\033[F\033[J";
       }

       if( orf_dictionary.find(match.query) == orf_dictionary.end() )  {
          orf_dictionary[match.query] = match_vector;
       }

       orf_dictionary[match.query].push_back(match);
       _num_orfs += 1;
 //      std::cout << match.subject << "  " << match.query << std::endl;

    }
   std::cout << "done" << std::endl;

   unsigned int maxReadLength = getMaxReadSize(orf_dictionary, contigs_dictionary); 

     int j = 0;

     std::map<string, vector<MATCH> >::iterator itcont;
     vector<MATCH>::iterator itorf;
     for(itcont= orf_dictionary.begin(); itcont != orf_dictionary.end(); itcont++) 
     {
         for(itorf=orf_dictionary[itcont->first].begin(); itorf != orf_dictionary[itcont->first].end(); itorf++)  {
          j++;
          if( j%1000==0) {
            
           std::cout << "\n\033[F\033[J" << j;
           //std::cout << "\n\033[F\033[J";
         }

     //  if( orfnames.find(match.subject) != orfnames.end() ) {
           try {
           //    std::cout << match.query << std::endl;
               substring_coverage(contigs_dictionary, itorf->query, itorf->start, itorf->end, coverage, maxReadLength); 
               orf_length += itorf->end - itorf->start;
//               std::cout << coverage.coverage << "  " << coverage.numreads << "  " << coverage.substring_length << "   " << coverage.uncovered_length << std::endl;

                if(count_type==0) {
                   orfnames[itorf->subject] = (1E9/static_cast<float>(num_mappable_reads))*(static_cast<float>(coverage.numreads)/static_cast<float>(coverage.substring_length));
                } else if(count_type==1) {
                   orfnames[itorf->subject] = static_cast<float>(coverage.numreads);
                } else if(count_type==2){
                   orfnames[itorf->subject] = (1E3/genome_equivalent)*(static_cast<float>(coverage.numreads)/static_cast<float>(coverage.substring_length));
                }
              // std::cout << num_mappable_reads << "  " << coverage.numreads << "  " <<  match.subject << "  " << orfnames[match.subject] << std::endl;
           }

           catch(...) {
               std::cout << "error\n";
               orfnames[match.subject] = 0;
               exit(0);
           }
      // }
       }

       //std::cout << match.query << " " << match.subject << "  " << match.start << " " << match.end << std::endl;

    }
   return _num_orfs;
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

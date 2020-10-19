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

void read_bam_files(CONTIG_ORF *contig_orf, vector<string> read_map_files) {

  SamFile samIn;
  SamFileHeader samFileHeader;
  SamRecord samRecord;

  for (size_t i = 0; i < read_map_files.size(); i++) {
     samIn.OpenForRead(read_map_files[i].c_str());
     samIn.ReadHeader(samFileHeader);
     while (samIn.ReadRecord(samFileHeader, samRecord)) {
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
}

void read_orf_names(string pathways_table_filename, map<string, float> &orfnames) {

    char buf[1000000];
    std::ifstream input(pathways_table_filename.c_str());
    if(!input.good()){
        std::cerr << "Error opening '"<<pathways_table_filename<<"'. Bailing out." << std::endl;
        return ;
    }   

    string stringCOMMENT("PWY_NAME");
    vector<char *> fields; 

    std::string line;
    while( std::getline( input, line ).good() ){
       if( matchString(line, stringCOMMENT, true) ) continue;
       fields.clear();
       split(line, fields, buf); 
       if( fields.size() <= 5) continue; 

       for(vector<char *>::iterator it=fields.begin()+5; it!=fields.end() ; it++) {
           orfnames[std::string(*it)] = 0;
       }
    }

    input.close();

}


RUN_STATS  detect_multireads_blastoutput(const std::string &blastoutput_file, const std::string &format,\
     vector<MATCH> &all_reads, map<std::string, unsigned long> &multireads, bool show_status) {

    MatchOutputParser *parser = ParserFactory::createParser(blastoutput_file, format);
    if( parser ==0 ) {
        std::cout << "ERROR : Cannot open a parser to parse the file " << blastoutput_file << std::endl;
    }

    MATCH  match;
    map<std::string, unsigned long> _multireads;

    map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned  int > > reads_dict;
   // map<std::string, std::pair<bool, bool > > reads_dict;
 //   map<std::string, std::pair<unsigned int, unsigned int > > multi_reads_dict;
   
    if( show_status) std::cout << std::endl << "Number of reads processed : " ;

    RUN_STATS stats;

    int i ;
    struct QUADRUPLE <bool, bool, unsigned int, unsigned int>   p;
    for( i =0; ; i++ )  {
       if( !parser->nextline(match) )  break;
       if( i >= _MAX ) break;

       if( match.mapped)  stats.num_mapped_reads++;
       else  stats.num_unmapped_reads++;

       if( match.parity)  stats.num_reads_2++; else stats.num_reads_1++;
       //if( match.multi)  num_multireads++;

       if( reads_dict.find(match.query) == reads_dict.end()) {
            p.first = false; 
            p.second = false;
            p.third = 0; 
            p.fourth = 0;
            reads_dict[match.query] = p;
       }
       stats.num_total_reads++;
      
       // if it is not mapped then ignore it
       if( !match.mapped)  continue;

       if( match.parity ) {
          reads_dict[match.query].first = true;
          reads_dict[match.query].third++;
       }
       else {
          reads_dict[match.query].second = true;
          reads_dict[match.query].fourth++;
       }

       // store it to process later by looking up the dictionary

       try {
          all_reads.push_back(match);
       }
       catch(...) {
          cout << "failing " << match.query << "   " << all_reads.size() <<  endl;
       }

       if( show_status && i%10000==0) {
           std::cout << "\n\033[F\033[J";
           std::cout << i ;
       }
       //std::cout << match.query << "   " << match.subject <<  " "  << match.start << " " << match.end << std::endl;
    }

    for( map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> >::iterator it = reads_dict.begin(); it != reads_dict.end(); it++)  {
        if( !(it->second.first && it->second.second))  stats.num_singleton_reads++;
        if( it->second.third > 1)  { stats.num_multireads++; stats.num_secondary_hits += it->second.third -1 ; }
        if( it->second.fourth  > 1) {  stats.num_multireads++; stats.num_secondary_hits += it->second.fourth-1; }
    }


    stats.num_distinct_reads_unmapped = stats.num_unmapped_reads;
    stats.num_distinct_reads_mapped = stats.num_mapped_reads - stats.num_secondary_hits;


    for( vector<MATCH>::iterator it = all_reads.begin(); it != all_reads.end(); it++)  {

       if( it->parity == 0  ) {
           if( reads_dict[it->query].first && reads_dict[it->query].second )
               it->w = 0.5/static_cast<float>(reads_dict[it->query].third);
           else
               it->w = 1/static_cast<float>(reads_dict[it->query].third);
       }
       else  { //parity 1
           if( reads_dict[it->query].first && reads_dict[it->query].second )
               it->w = 0.5/static_cast<float>(reads_dict[it->query].fourth);
           else
               it->w = 1/static_cast<float>(reads_dict[it->query].fourth);
       }
//    all_sequence_reads += n;
    }
    // now store the multireads into the multireads map variable
    
    //*_num_unmapped_reads= parser->get_Num_Unmapped_Reads() ;
    delete parser;
    return stats;
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



void add_RPKM_value_to_pathway_table(const string &pathways_table_filename, const string &output_file, map<string, float> &orfnames) {
    char buf[1000000];
    std::ifstream input(pathways_table_filename.c_str());
    if(!input.good()){
        std::cerr << "Error opening '"<<pathways_table_filename<<"'. Bailing out." << std::endl;
        return ;
    }   

    string stringCOMMENT("PWY_NAME");

/*
    vector<char *> fields; 

       fields.clear();
       split(line, fields, buf); 
       if( fields.size() <= 5) continue; 
*/
    vector<string> lines;
    std::string line;

    //read in all the lines
    string headerline;
    while( std::getline( input, line ).good() ){
       if( matchString(line, stringCOMMENT, true) ) { headerline = line;  continue;}
       lines.push_back(line);
    }
    input.close();

    std::ostream *output  = &std::cout;
    std::ofstream fout;
    if( output_file.size()!=0) {
        fout.open(output_file.c_str(), std::ifstream::out);
        output = &fout;
    }

    
    float pwy_rpkm;
    vector<char *> fields; 
    *output << headerline << "\tRPKM_COUNT" << std::endl ;
    for(vector<string>::iterator it=lines.begin(); it != lines.end(); it++) {
        split(*it, fields, buf, '\t'); 
        if( fields.size() < 5) continue; 
         
        for(int i =0; i< 5; i++)  {
            if( i > 0 ) *output << '\t';  
            *output << fields[i]  ;
        }

       
       pwy_rpkm = 0;
       for(vector<char *>::iterator it=fields.begin()+5; it!=fields.end() ; it++) {
         // std::cout << *it << std::endl;
          pwy_rpkm += orfnames[std::string(*it)];
       }
       sprintf(buf, "\t%0.2f", pwy_rpkm) ;
       *output << buf;

       for(vector<char *>::iterator it=fields.begin()+5; it!=fields.end() ; it++) {
          *output << "\t" << *it;
       }
       *output << std::endl ;
    }
    if( output_file.size() !=0 )  fout.close();
    
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

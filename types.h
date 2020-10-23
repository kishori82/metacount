#ifndef __DATA_TYPES__
#define __DATA_TYPES__

#include <map>
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

typedef struct _ORFINFO {
  std::string id;
  uint32_t start, end;
  float count;
  
  _ORFINFO(string _id, uint32_t _start, 
           uint32_t _end): id(_id), start(_start), 
           end(_end), count(0) {}
} ORFINFO;

enum IDTYPE {ORFID, CONTIGID};

typedef map<std::string, vector<ORFINFO *> *>  CONTIG_ORF;

typedef struct _TRIPLET {
    unsigned int start, end;
    float multi;
} TRIPLET;

typedef vector<TRIPLET> TRIPLETS;

typedef struct _CONTIG {
    unsigned int L;
    TRIPLETS M;
} CONTIG;

enum READNO {FIRST, SECOND};
enum ALNVALID {VALID, INVALID};
enum READTYPE {SE, PE};
enum ALIGNTYPE {NO_ALIGN, PRIMARY, NON_PRIMARY, OTHER};
enum CHIMERIC{SUPPLEMENTARY,  NON_SUPPLEMENTARY};

typedef  struct _MATCH {
    std::string query, subject;
    int16_t start, end, alnlength;

    int16_t samflag;

    enum ALNVALID alnvalid;   // VALID, INVALID
    enum READNO readno;   //  FIRST, SECOND
    enum READTYPE readtype;   //;  SE, PE (single-end or paired-end)    
    enum ALIGNTYPE aligntype;  // PRIMARY, NON_PRIMARY
    enum CHIMERIC chimeric;  // SUPPLEMENTARY, NON_SUPPLEMENTARY
 
    float weight;

    _MATCH(): weight(0) { } 

    string print() {
        string out(this->query);
        out = out + "\t";
        out = out + std::to_string(this->samflag);
        out = out + "\t";
        out = out + this->subject;
        out = out + "\t";
        out = out + std::to_string(this->start);
        out = out + "\t";
        out = out + std::to_string(this->end);
        out = out + "\t";
        out = out + std::to_string(this->alnlength);
        out = out + "\t";
        if (this->readtype == SE) {
            out = out + "SE";
        } else { 
            out = out + "PE";
        }
        out = out + "\t";
        if (this->readno == FIRST) {
            out = out + "R1";
        } else { 
            out = out + "R2";
        }
        out = out + "\t";
        if (this->alnvalid == VALID) {
            out = out + "VALID";
        } else { 
            out = out + "INVALID";
        }
        out = out + "\t";
        switch (this->aligntype) {
            case NO_ALIGN:
                out = out + "\t" + "NO_ALIGN";
                break;
            case PRIMARY:
                out = out + "\t" + "PRIMARY";
                break;
            case NON_PRIMARY:
                out = out + "\t" + "NON_PRIMARY";
                break;
            default:
                out = out + "\t" + "OTHER";
                break;
        } 

        switch (this->chimeric) {
            case SUPPLEMENTARY:
                out = out + "\t" + "SUPPLEMENTARY";
                break;
            case NON_SUPPLEMENTARY:
                out = out + "\t" + "NON_SUPPLEMENTARY";
                break;
            default:
                out = out + "\t" + "OTHER";
        }

        out = out + "\t" +  std::to_string(this->weight);
        return out;
    } 
} MATCH;


typedef struct _GLOBAL_STATS {
    int num_file; 
    int num_contigs; 
    int total_contig_lengths; 
    int total_reads1; 
    int total_reads2; 
    float tot_count_across_orf;

    _GLOBAL_STATS() {
      num_file = 0; 
      num_contigs = 0; 
      total_contig_lengths = 0; 
      total_reads1 = 0; 
      total_reads2 = 0; 
      tot_count_across_orf = 0;
    }
    void print_stats(std::ostream *output) {
        *output << std::endl;
        *output << "Number of sam/bam files           : " << num_file << std::endl; 
        *output << "Number of contigs                 : " << num_contigs << std::endl; 
        *output << "Total contig length               : " << total_contig_lengths << std::endl; 
        *output << "Total read1s                      : " << total_reads1 << std::endl;
        *output << "Total read2s                      : " << total_reads1 << std::endl;
        *output << "Total count to across orfs        : " << tot_count_across_orf << std::endl; 
    }
} GLOBAL_STATS;

typedef struct _RUN_STATS {
    int num_unmapped_alns ; 
    int num_mapped_alns ; 
    int num_reads_1 ; 
    int num_reads_2 ; 
    int num_reads1_on_orfs ; 
    int num_reads2_on_orfs ; 
    int num_total_alns ;
    
    _RUN_STATS() {
        num_unmapped_alns =0; 
        num_mapped_alns =0; 
        num_reads_1 =0; 
        num_reads_2 =0; 
        num_total_alns = 0;
        num_reads1_on_orfs = 0 ; 
        num_reads2_on_orfs = 0; 
    }

    struct _RUN_STATS&  operator+( const struct _RUN_STATS & stats) {
        this->num_unmapped_alns += stats.num_unmapped_alns;  
        this->num_mapped_alns  += stats.num_mapped_alns  ;
        num_reads_1 += stats.num_reads_1;  
        num_reads_2 += stats.num_reads_2; 
        num_total_alns  += stats.num_total_alns  ;
        return *this;
    }


    void print_stats(std::ostream *output) {
        *output << std::endl;
        *output << "Number of alignments             : " << num_total_alns << std::endl; 

        int x = num_total_alns == 0 ? 1 : num_total_alns;           

        *output << "Number of unmapped alignments    : " << num_unmapped_alns << std::endl;
        *output << "Number of mapped alignments      : " << num_mapped_alns << std::endl; 

        *output << "Percentage of mapped alignments  : " << std::setprecision(4)  
                << static_cast<float>(num_mapped_alns) *100/static_cast<float>(x) << "%" <<  std::endl; 
        *output << "Number of read 1                 : " << num_reads_1 << std::endl; 
        *output << "Number of read 2                 : " << num_reads_2 << std::endl; 
        *output << "Number of read1 on ORFs          : " << num_reads1_on_orfs << std::endl; 
        *output << "Number of read2 on ORFs          : " << num_reads2_on_orfs << std::endl; 
    }
} RUN_STATS;

typedef struct _RESULTS {
   vector<RUN_STATS *> runresults;
   GLOBAL_STATS global_stats;
   CONTIG_ORF *contig_orf;
   std::map<string, uint32_t> *contig_lengths;
} RESULTS;

#endif //__DATA_TYPES__




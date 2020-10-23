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
  
  _ORFINFO(string _id, uint32_t _start, uint32_t _end): id(_id), start(_start), end(_end), count(0) {}
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

template< typename A, typename B, typename C, typename D>
struct QUADRUPLE {
     A first;
     B second;
     C third;
     D fourth;
};

typedef struct _RUN_STATS {
    int num_unmapped_alns ; 
    int num_mapped_alns ; 
    int num_singleton_reads ; 
    int num_reads_1 ; 
    int num_reads_2 ; 
    int num_multireads ;
    int num_total_alns ;
    int num_secondary_hits ;
    int num_distinct_reads_unmapped;
    int num_distinct_reads_mapped;
    
    _RUN_STATS() {
        num_unmapped_alns =0; 
        num_mapped_alns =0; 
        num_singleton_reads =0; 
        num_reads_1 =0; 
        num_reads_2 =0; 
        num_multireads = 0;
        num_total_alns = 0;
        num_secondary_hits = 0;
        num_distinct_reads_unmapped =0 ;
        num_distinct_reads_mapped = 0;
    }

    struct _RUN_STATS&  operator+( const struct _RUN_STATS & stats) {
        this->num_unmapped_alns += stats.num_unmapped_alns;  
        this->num_mapped_alns  += stats.num_mapped_alns  ;
        num_singleton_reads += stats.num_singleton_reads ;
        num_reads_1 += stats.num_reads_1;  
        num_reads_2 += stats.num_reads_2; 
        num_multireads += stats.num_multireads; 
        num_total_alns  += stats.num_total_alns  ;
        num_secondary_hits  += stats.num_secondary_hits;
        num_distinct_reads_unmapped  += stats.num_distinct_reads_unmapped  ;
        num_distinct_reads_mapped  += stats.num_distinct_reads_mapped  ;
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
        //*output << "Number of singletons             : " << num_singleton_reads << std::endl; 
        *output << "Number of read 1                 : " << num_reads_1 << std::endl; 
        *output << "Number of read 2                 : " << num_reads_2 << std::endl; 
        //*output << "Number of multireads             : " << num_multireads << std::endl; 
        //*output << "Number of secondary hits         : " << num_secondary_hits << std::endl; 
       // *output << "Number distinct reads mapped     : " << num_distinct_reads_mapped << std::endl; 
        //*output << "Number distinct reads unmapped   : " << num_distinct_reads_unmapped << std::endl; 
    }
} RUN_STATS;

#endif //__DATA_TYPES__




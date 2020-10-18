#ifndef __UTILITIES__H__
#define __UTILITIES__H__

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <cstdlib>
#include <regex>

#include "types.h"

using namespace std;

// Structure for RPKM input options
struct Options {
    /* Input files for RPKM */
    string contigs_file; // the contigs file
    string orf_file; // .gff ORF file from PRODIGAL
    string pathways_table; // a table from Pathway Tools with ORFs
    string output_file; // location to write output file (i.e., update pathway table)
    string stats_file; // location to write output file (i.e., update pathway table)
    string orf_rpkm_file; // location to write output file (i.e., update pathway table)
    vector<string> read_map_files;
    
    /* Flags and settings */
    bool multi_reads; // flag for detecting multiple mapping of reads
    bool show_status; // shows the counter that countes the number of reads processed, and other info 
                       // on screen
    bool read_counts; // show counts not rpkm
    string reads_map_file_format; // aligner type BWA or BLAST, two SAM files or one
    
    int count_type;
    float genome_equivalent;   //genome equivanet from RPKG such as MicrobeSensus
    // Constructor with default settings 
    Options(){
       contigs_file = "";
       stats_file = ""; // location to write output file (i.e., update pathway table)
       read_map_files.clear();
       orf_file = "";
       pathways_table = "";
       output_file = "";
       output_file = "";
       orf_rpkm_file = ""; // location to write output file (i.e., update pathway table)

       multi_reads = false;
       show_status = false;
       read_counts = false;
       reads_map_file_format = "sam-1";  
       genome_equivalent = 1;
       count_type = 0;
    };
    
    void print_usage( char *arg);
    void print_options();
 
    bool SetOptions(int argc, char *argv[]);
	// bool SetOptionCommon( const char *flag, const char *value );
	// bool SetOption( const char *flag, const char *value );
	// bool SetOption2D( const char *flag, const char *value );
	// bool SetOptionEST( const char *flag, const char *value );
	// bool SetOptions( int argc, char *argv[], bool twodata=false, bool est=false );

	void Validate();
	// void ComputeTableLimits( int naan, int typical_len, size_t mem_need );

	void Print();
};

void split(const std::string  &strn, std::vector<char *> &v, char *buf, char d='\t');
std::string get_orf_name(std::string & strn, std::vector<char *> &v, char *buf);


bool matchString(const string &str, const string & stringtomatch, bool fromstart=false);

void get_fasta_sequence_info(const std::string &fasta_file_name);

std::string extract_sequence_name(const std::string &name);

string to_string(unsigned long i);


// Print system error and exit
void error(char *msg);

string shorten_id(string name, enum IDTYPE);

uint16_t strToUint16(string str1);

#endif // __UTILITIES__H__


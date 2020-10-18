#ifndef __HELPER_H__
#define __HELPER_H__
#include <map>
#include <string>
#include <ostream>
#include <iterator>
#include <assert.h>
#include "types.h"
#include "parser.h"

using namespace std;
CONTIG_ORF * create_contig_orf_map(const string &orf_file);

void read_orf_names(std::string pathways_table, std::map<string, float> &orfnames) ;

unsigned long create_contigs_dictionary(std::string contigs_file, 
                                  std::map<std::string, CONTIG> &contigs_dictionary);

RUN_STATS  detect_multireads_blastoutput(const std::string &blastoutput_file,\
                             const std::string &format, vector<MATCH> &all_reads,\
                             map<std::string, unsigned long> &multireads,\
                              bool show_progress_counter = false);

void  process_blastoutput(const std::string & reads_map_file, std::map<string,\
                           CONTIG> &contigs_dictionary,\
                           const std::string &reads_map_file_format,\
                            vector<MATCH> &all_reads, map<string,\
                            unsigned long> & multireads, bool show_status= false);

void substring_coverage(std::map<string, CONTIG> &contigs_dictionary, const std::string &contig,\
                         unsigned long start, unsigned long end, COVERAGE &coverage, unsigned int maxReadLength);

unsigned long ORFWise_coverage( map<string, CONTIG> &contigs_dictionary, const string &orf_file,\
                              map<string, float> &orfnames, unsigned long genome_length,\
                              unsigned long &orf_length,  unsigned long num_mappable_reads,\
                              int count_type, float genome_equivalent= 1, bool show_status = false);


void add_RPKM_value_to_pathway_table(const string &pathways_table, const string &output_file, map<string, float> &orfnames);

void writeOut_ORFwise_RPKM_values(const string orf_rpkm_file,  map<string, float> &orfnames);

#endif // __HELPER_H__

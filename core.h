#ifndef __HELPER_H__
#define __HELPER_H__
#include <map>
#include <string>
#include <ostream>
#include <iterator>
#include <assert.h>
#include "types.h"

#include "parser.h"
#include "SamValidation.h"

using namespace std;
CONTIG_ORF * create_contig_orf_map(const string &orf_file);

std::map<string, uint32_t> * read_bam_files(CONTIG_ORF *contig_orf, vector<string> read_map_files);

std::map<string, uint32_t> * get_contig_information(vector<string> read_map_files);

unsigned long create_contigs_dictionary(std::string contigs_file, 
                                  std::map<std::string, CONTIG> &contigs_dictionary);

void writeOut_ORFwise_RPKM_values(const string orf_rpkm_file,  map<string, float> &orfnames);

#endif // __HELPER_H__

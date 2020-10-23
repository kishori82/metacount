#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <experimental/filesystem>
#include <functional>
#include <numeric>


#include "types.h"
#include "utilities.h"
#include "options.h"
#include "core.h"
#include "SamValidation.h"

using namespace std;
namespace fs = std::experimental::filesystem;


bool compare_triplets(const TRIPLET &a, const TRIPLET &b) {
  return a.start < b.start ? true : false; 
}
 
int main(int argc, char **argv ){
  // parse options
  Options options;

  // if (argc < 9) { options.print_usage(argv[0]); exit(1); }
  if (options.SetOptions(argc, argv)==false) { 
    options.print_usage(argv[0]); exit(0); 
  }

  RESULTS results;

  results.contig_orf = create_contig_orf_map(options.orf_file);

  //CONTIG_ORF *contig_orf = create_contig_orf_map(options.orf_file);

  for (size_t i = 0; i < options.read_map_files.size(); i++) {
     if (!fs::exists(options.read_map_files[i].c_str())) {
        std::cout << "ERROR: file" << options.read_map_files[i] << " is missing!\n";
        std::cerr << "ERROR: File " << options.read_map_files[i] << " is missing!\n";
        exit(0);
     }
  }

  std::map<string, uint32_t> *contig_lengths = get_contig_information(options.read_map_files);

  std::cout << "num contigs " << contig_lengths->size() << std::endl;

  //print_contig_orf_map(contig_orf);

  sort_contig_orf_map(results.contig_orf);

  std::map<string, uint32_t> * contig_read_counts = 
    read_bam_files(results.contig_orf, options.read_map_files);

  //print_contig_read_counts(contig_read_counts);

  //print_contig_orf_map(contig_orf);

  uint32_t genome_length = accumulate(
                                      contig_lengths->begin(), 
                                      contig_lengths->end(), 0, 

                                      [] (uint32_t value, const std::pair<string, uint32_t>&p) {
                                          return value + p.second; 
                                      }
                                     );

                 //[] (uint32_t value, const std::pair<string, vector<ORFINFO *>*> &orfvec) {
  uint32_t total_read_count_across_orfs
    = accumulate(
                 results.contig_orf->begin(), 
                 results.contig_orf->end(), 0, 
                    
                 [] (uint32_t value, const std::pair<string, vector<ORFINFO *>*> &orfvec) {
                       float tot_count = value;
                       for (auto it = orfvec.second->begin(); it != orfvec.second->end(); it++) {
                          tot_count = tot_count + (*it)->count; 
                       }
                       return tot_count; 
                 }
               );




  std::cout << "Total genome length " << genome_length << std::endl;
  std::cout << "Total read count across all orfs " << total_read_count_across_orfs << std::endl;
  exit(0);
   
  //options.print_options();

//  std::cout << "number of contigs " << contigs_dictionary.size() << std::endl;

  std::ostream *output;
  std::ofstream fout;
  if (options.stats_file.size() > 0) {
        fout.open(options.stats_file.c_str(), std::ifstream::out);
        output = &fout;
    }   
 
  if (options.show_status) { 
       std::cout << "Composite stats for all files " << std::endl;
  }

  if( options.show_status ) std::cout << "\n\nSorting  the read matches .....";



  return 0;
} 

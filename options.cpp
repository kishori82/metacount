#include "options.h"
//#include <iostream>
//#include <fstream>


void Options::print_usage(char *arg) {
  std::cout << "USAGE : "   << arg\
             << "      : -h [ for help ]\n"\
             << "      : --sam  <sam/bam> \n"\
             << "      : --gff  <orf_file_gff> [REQUIRED]\n"
             << "      : --stats-output-file stats_file [OPTIONAL]\n"\
             << "      : --o <outputfile> \n"\
             << "      : --stats [prints overall stats]\n"\
             << "      : --stats-type [TPM/RPKM/COUNT/ALL default TPM]\n"\
             << std::endl;
}

bool Options::SetOptions(int argc, char *argv[]) { 
    for(int i = 1; i < argc ; i++) {   
        if (strncmp(argv[i], "-h", strlen("-h")) == 0) {
            print_usage(argv[0]);
            exit(0);
        }   
        else if (strncmp(argv[i], "--stats-type", strlen("--stats-type")) == 0) {
            this->count_type = argv[++i];
            if (this->count_type != "ALL" && this->count_type != "TPM" 
                && this->count_type != "RPKM" && this->count_type != "COUNT") {
               print_usage(argv[0]);
               exit(0);
            }
        }   
        else if (strncmp(argv[i], "--stats-output-file", strlen("--stats-output-file")) == 0) {
            this->stats_file = argv[++i];
        }   
        else if (strncmp(argv[i], "--sam", strlen("--sam")) == 0) {
            this->read_map_files.push_back(string(argv[++i]));
        }   
        else if (strncmp(argv[i], "--stats", strlen("--stats")) == 0) {
            this->show_stats = true;
        }   
        else if (strncmp(argv[i], "--o", strlen("--o")) == 0) {
            this->output_file = argv[++i];
        } 
        else if (strncmp(argv[i], "--gff", strlen("--gff")) == 0) {
            this->orf_file = argv[++i];
        }
        else {
            std::cout << "ERROR: Input optiono" << std::endl;
            print_usage(argv[0]);
            exit(0);
        }
    } //for loop for arguments processing

    if (this->read_map_files.size() == 0) {
        std::cout << "ERROR: There must be a at least one read map file" << std::endl;
        return false;
    }
    
    if (this->orf_file.size()==0) {
        std::cout << "ERROR: There must be an ORF file [ GFF format ]" << std::endl;
        return false;
    }

    return true;
};

void Options::print_options() {
    std::cout << "Alignment file (SAM format)  : "<< this->read_map_files.size() << std::endl; 
    std::cout << "ORF file [GFF]                    : "<< this->orf_file << std::endl; 
    std::cout << "Output file                       : "<< this->output_file << std::endl; 
}


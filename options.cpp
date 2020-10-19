#include "options.h"
//#include <iostream>
//#include <fstream>


void Options::print_usage(char *arg) {
  std::cout << "USAGE : "   << arg\
             << "      : -h [ for help ]\n"\
             << "      : --sam  <sam/bam> \n"\
             << "      : --gff  <orf_file_gff> [REQUIRED]\n"
             << "      : --ORF-RPKM  <orf_rpkm_file> [default: None]\n"\
             << "      : --stats stats_file [OPTIONAL]\n"\
             << "      : --o <outputfile> \n"\
             << "      : --m <multireads>  [OPTIONAL , default on]\n"\
             << "      : --status [shows status]\n"\
             << "      : --type [rpkm[0]/count[1]/rpkg[2] default 0 ]\n"\
             << "      : --geonme_equivalent [genome equivalent from rpkg]\n"
             << "      : --count  [shows counts]"\
             << std::endl;
}

bool Options::SetOptions(int argc, char *argv[]) { 
    for(int i = 1; i < argc ; i++) {   
        if (strncmp(argv[i], "-h", strlen("-h")) == 0) {
            print_usage(argv[0]);
            exit(0);
        }   
        else if (strncmp(argv[i], "--type", strlen("--type")) == 0) {
            this->count_type = atoi(argv[++i]);
        }   
        else if (strncmp(argv[i], "--genome_equivalent", strlen("--genome_equivalent")) == 0) {
            this->genome_equivalent =  atof(argv[++i]);
        }   
        else if (strncmp(argv[i], "--stats", strlen("--stats")) == 0) {
            this->stats_file = argv[++i];
        }   
        else if (strncmp(argv[i], "--sam", strlen("--sam")) == 0) {
            this->read_map_files.push_back(string(argv[++i]));
        }   
        else if (strncmp(argv[i], "--status", strlen("--status")) == 0) {
            this->show_status = true;
        }   
        else if (strncmp(argv[i], "--o", strlen("--o")) == 0) {
            this->output_file = argv[++i];
        } 
        else if (strncmp(argv[i], "--gff", strlen("--gff")) == 0) {
            this->orf_file = argv[++i];
        }
        else if (strncmp(argv[i], "--ORF-RPKM", strlen("--ORF-RPKM")) == 0) {
            this->orf_rpkm_file = argv[++i];
        }
        else if (strncmp(argv[i], "--m", strlen("--m")) == 0) {
            this->multi_reads = true;
        }
        else {
            std::cout << "ERROR: Choices for -format option must of type blastout, sam-1 or sam-2" << std::endl;
               return false;
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
    std::cout << "Multi read treatment              : "<< this->multi_reads << std::endl; 
}


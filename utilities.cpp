#include "utilities.h"
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


std::string extract_sequence_name(const std::string &name) {
     char  cstr[1000000];
     strcpy(cstr, name.c_str());
     
     char * cptr = cstr;

     while( *cptr != '\t' && *cptr !=  ' ' && *cptr != '\0' )  cptr++; 
     (*cptr) ='\0';

     std::string sname(cstr);
     return sname;
}



void split(const string  &strn, std::vector<char *> &v, char *buf, char d) {
  strcpy(buf, strn.c_str());
  char *s1 = buf;
  v.clear();
  v.push_back(s1);
  while(*s1 != '\0') {
     if(*s1==d) { 
       *s1 = '\0';
       v.push_back(s1+1);
     }
     s1++;
  }
}

std::string get_orf_name(std::string  &strn, std::vector<char *> &v, char *buf) {
    split(strn, v, buf, ';'); 
    if(v.size() == 0)  return std::string("");
    split(std::string(v[0]), v, buf, '=');
    if(v.size() < 2)  return std::string("");
    return std::string(v[1]);
}

bool matchString(const string &str, const string & stringtomatch, bool fromstart) {
    unsigned long pos = str.find(stringtomatch);
    if(fromstart &&  pos ==0 ) return true;

    if( !fromstart && pos >= 0) return true;
    return false;

}

// Print system error and exit
void error(char *msg) {
    perror(msg);
    exit(1);
}

string shorten_id(string name, enum IDTYPE idtype)  {
  regex contigid_regex("_([0-9]+)$");
  regex orfid_regex("_([0-9]+_[0-9]+)$");
  std::smatch sm; 

  switch (idtype) {
     case ORFID: 
       std::regex_search(name, sm, orfid_regex);
       break;
     case CONTIGID: 
       std::regex_search(name, sm, contigid_regex);
       break;
     default: 
       break;
  } 

  if (sm.size() > 1) {
     return sm[1];
  }

  return string("");
}

uint16_t strToUint16(string str1) {

  int myint(std::stoi(str1));

  uint16_t myint16(0);
  if (myint <= static_cast<int>(UINT16_MAX) && myint >=0) {
     myint16 = static_cast<uint16_t>(myint);
  } else {
     perror("ERROR : Cannot convert orf start position to 16 bit  integer\n");
  }
  return myint16;
}


#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include "genes.h"

void Genes::readFastas(std::string directory) {

    for (const auto & entry : std::filesystem::directory_iterator(directory)) {
        if(entry.is_directory()) {
          std::cerr << "\33[2K\rreading genes of " << entry.path().filename().string() << std::flush; // white space to write over longer names
          std::string fastafile = "";
          for (const auto & file : std::filesystem::directory_iterator(entry)) {
              if(file.is_regular_file() && file.path().extension() == ".fasta") {
                  // std::cerr << "fasta file: " << file.path() << std::endl;
                  // readFasta("bdi", "/data/bls/newdata/nextflow_output/bdi/masked_promoter_bdi.fasta");
                  if (fastafile.empty()) fastafile = file.path().string();
                  else {
                      if(file.path().string().find("geneid") != std::string::npos)
                        fastafile = file.path().string();
                  }
              }
          }
          readFasta(entry.path().filename(), fastafile);
        }
        // std::cerr << "genemap size: " << genemap.size() << std::endl;
    }
    // for (auto x : genemap) {
        // std::cerr << x.first << " -> " << x.second << std::endl;
    // }
    std::cerr << '\r' << genemap.size() << " genes in genemap" << std::endl;
}
Gene *Genes::getGene(std::string geneid) {
    if(genemap.find(geneid) != genemap.end())
        return &genemap[geneid];
    else
        std::cerr << "geneid " << geneid << " not found in genemap" << std::endl;
    return NULL;
}
char Genes::getRandomCharFromIupac(char c) {
    char new_c = characterToMask[c][rand() % characterToMask[c].size()];
    // std::cerr << "replaced " << c << " with " << new_c << std::endl;
    return new_c;
}
void Genes::fixSequence(std::string &seq) {
    std::for_each(seq.begin(), seq.end(), [](char & c) { // convert all to upper case!
        c = ::toupper(c);
        if (validCharacters.find(c) == validCharacters.end()) {
            c = getRandomCharFromIupac(c);
        }
    });
}


const std::unordered_set<char> Genes::validCharacters ({ 'A', 'C', 'G', 'T', 'N'  });
const std::vector<std::string> Genes::characterToMask ({
    "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", // 0
    "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", // 16
    "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", // 32
    "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", // 48
    "", "A", "CGT", "C", "AGT", "", "", "G", "ACT", "", "", "GT", "", "AC", "ACGT", "", // 64
    "", "", "AG", "GC", "T", "", "ACG", "AT", "", "CT" // 80
});
void Genes::readFasta(std::string species, std::string fasta) {
    // read contents of fasta and add to the map

    // std::cerr<<" reading species: " << species << " from " << fasta << std::endl;
    std::ifstream f(fasta);
    std::string id, parsedid, seq;
    if (f.is_open()) {
        while ( getline (f, id) )  {
            if(!id.empty() && id[0] == '>' ) {
                parsedid = id.substr(1, id.find_first_of(':') - 1);
                if(getline(f, seq)) {
                    if(!seq.empty()) {
                        fixSequence(seq);
                        genemap[parsedid] = {species, seq};
                    } else {
                        std::cerr << "empty seq line..." << std::endl;
                    }
                } else {
                    std::cerr << "incorrect fasta format, no sequence line after id line" << std::endl;
                }
            } else {
                std::cerr << "incorrect fasta format, no id line found" << std::endl;
            }
        }
        f.close();
    }


}


#include <string>
#include <vector>
#include <fstream>
#include "genefamily.h"

std::chrono::time_point<std::chrono::system_clock> prevTime;

const std::unordered_set<char> GeneFamily::validCharacters ({ 'A', 'C', 'G', 'T', 'N', ' ', '$'  });

void startChrono()
{
        prevTime = std::chrono::system_clock::now();
}

double stopChrono()
{
        std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - prevTime;
        return (elapsed.count());
}

void GeneFamily::readOrthologousFamily(const int mode, const std::string& filename, const std::vector<float> blsThresholds_, const Alphabet alphabet,
const int type, const std::pair<short, short> l, const int maxDegeneration) {
    std::ifstream ifs(filename.c_str());
    readOrthologousFamily(mode, ifs, blsThresholds_, alphabet, type, l, maxDegeneration);
}


void GeneFamily::readOrthologousFamily(const int mode, std::istream& ifs, const std::vector<float> blsThresholds_, const Alphabet alphabet,
const int type, const std::pair<short, short> l, const int maxDegeneration) {
  int totalCount = 0;
  while (ifs) {
    std::vector<size_t> stringStartPositions;
    std::vector<size_t> next_gene_locations;
    std::vector<std::string> gene_names;
    stringStartPositions.push_back(0);
    // READ DATA
    std::string T, newick, line, name;
    int N;
    getline(ifs, line);
    while(ifs && line.empty()) {getline(ifs, line);}
    if(!ifs || line.empty()) {continue;}
    name = line;
    getline(ifs, newick);
    getline(ifs, line);
    N = std::stoi(line);
    size_t current_pos = 0;
    next_gene_locations.push_back(current_pos);
    for (int i = 0; i < N; i++) {
        getline(ifs, line);
        // gene names
        std::vector<std::string> genes;
        line = line.substr(0, line.find_first_of('\t'));
        size_t start = 0;
        size_t end = line.find_first_of(' ', start);
        while(end != std::string::npos) {
            genes.push_back(line.substr(start, end - start));
            start = end + 1;
            end = line.find_first_of(' ', start);
        }
        genes.push_back(line.substr(start));
        for (size_t k =0; k < genes.size(); k++) {
            gene_names.push_back(genes[k]);
        } // add RC genes
        for (size_t k = genes.size() ; k > 0; k--) {
            gene_names.push_back(genes[k-1]);
        }
        // genes
        getline(ifs, line);
        if (!T.empty())
            T.push_back(IupacMask::DELIMITER);
        std::for_each(line.begin(), line.end(), [](char & c) { // convert all to upper case!
            if(c == IupacMask::FILLER)
                c = IupacMask::DELIMITER;
            else {
                c = ::toupper(c);
                // should already be fixed in the preprocessing tool!!!!
                // if (validCharacters.find(c) == validCharacters.end()) {
                    // c = IupacMask::getRandomChar(c);
                // }
            }
        });
        T.append(line);
        T.push_back(IupacMask::DELIMITER);
        stringStartPositions.push_back(T.size());
        T.append(Motif::ReverseComplement(line));
        stringStartPositions.push_back(T.size() + 1);
        // std::cout << T << std::endl;

        // add gene start locations...
        std::vector<size_t> gene_sizes;
        start = 0;
        end = line.find_first_of(' ', start);
        while(end != std::string::npos) {
            gene_sizes.push_back(end + 1 - start);
            start = end + 1;
            end = line.find_first_of(' ', start);
        }
        gene_sizes.push_back(line.size() + 1 - start);
        for (size_t k =0; k < gene_sizes.size(); k++) {
            current_pos += gene_sizes[k];
            next_gene_locations.push_back(current_pos);
        } // add RC genes
        for (size_t k = gene_sizes.size(); k > 0; k--) {
            current_pos += gene_sizes[k - 1];
            next_gene_locations.push_back(current_pos);
        }
    }
    T.push_back(IupacMask::DELIMITER);

    std::cerr << "T: " << T.length() << std::endl;
    std::cerr << "[" << name << "] " << N << " gene families" << std::endl;
    // PROCESS DATA
    startChrono();
    BLSScore bls(blsThresholds_, newick, N);
    // std::cout << T << std::flush;
    SuffixTree ST(T, true, stringStartPositions, gene_names, next_gene_locations);

    if (mode == 1) {
        int count = ST.matchIupacPatterns(ifs, std::cout, bls, maxDegeneration, l.second);
        double elapsed = stopChrono();
        std::cerr << "[" << name << "] " << count <<  " motifs located in " << elapsed << "s" << std::endl;
    } else if (mode == 0) {
        int count = ST.printMotifs(l, alphabet, maxDegeneration, bls, std::cout, type == 0); // 0 == AB, 1 is AF
        size_t iteratorcount = ST.getMotifsIteratedCount();

        totalCount += count;
        double elapsed = stopChrono();
        std::cerr << "\33[2K\r[" << name << "] iterated over " << iteratorcount << " motifs" << std::endl;
        std::cerr << "[" << name << "] counted " << count << " valid motifs in " << elapsed << "s" << std::endl;
    } else {
        std::cerr << "wrong mode given: " << mode << std::endl;
    }
  }
  if (mode == 0) std::cerr << "total motifs counted: " << totalCount << std::endl;
}

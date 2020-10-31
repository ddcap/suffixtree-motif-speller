
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

void GeneFamily::readOrthologousFamily(const int mode, const std::string& filename, std::vector<float> blsThresholds_, Alphabet alphabet, int type, std::pair<short, short> l, int maxDegeneration) {
    std::ifstream ifs(filename.c_str());
    readOrthologousFamily(mode, ifs, blsThresholds_, alphabet, type, l, maxDegeneration);
}

void GeneFamily::readOrthologousFamily(const int mode, std::istream& ifs, std::vector<float> blsThresholds_, Alphabet alphabet, int type, std::pair<short, short> l, int maxDegeneration) {
  int totalCount = 0;
  while (ifs) {
    std::vector<size_t> stringStartPositions;
    stringStartPositions.push_back(0);
    // READ DATA
    std::string T, newick, line, name;
    int N;
    getline(ifs, line);
    while(ifs && line.empty()) {getline(ifs, line);}
    if(!ifs || line.empty()) {continue;}
    name = line;
    std::cerr << "test [" << name << "] " << std::endl;
    getline(ifs, newick);
    getline(ifs, line);
    N = std::stoi(line);
    for (int i = 0; i < N; i++) {
        getline(ifs, line);
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
    }
    T.push_back(IupacMask::DELIMITER);

    std::cerr << "T: " << T.length() << std::endl;
    std::cerr << "[" << name << "] " << N << " gene families" << std::endl;
    // PROCESS DATA
    startChrono();
    BLSScore bls(blsThresholds_, newick, N);
    // std::cout << T << std::flush;
    SuffixTree ST(T, true, stringStartPositions);

    if (mode == 1) {
        getline(ifs, line);
        // collect in a global list and sort it after and then get family/genes easy!
        while (!line.empty() && line.compare("-") != 0 ) { // loop over motifs until empty line or line with - signaling the end
            occurence_bits occurence;

            std::vector<std::pair<int, int>> foundpos = (type == 0) ?
                 ST.matchIupacPatternWithPositions(line, bls, maxDegeneration, occurence) :
                 ST.matchIupacPattern(line, bls, maxDegeneration, occurence);

            std::cerr << "motif " << line << " with occurence: " << std::bitset<8>(occurence) << " and BLS: " << bls.getBLSScore(occurence) << std::endl;
            ST.printMotifPositions(std::cout, line, foundpos, line.length());
            getline(ifs, line);
        }
    } else if (mode == 0) {
        int count = ST.printMotifs(l, alphabet, maxDegeneration, bls, std::cout, type == 0); // 0 == AB, 1 is AF
        size_t iteratorcount = ST.getMotifsIteratedCount();

        totalCount += count;
        double elapsed = stopChrono();
        std::cerr << "\33[2K\r[" << name << "] iterated over " << iteratorcount << " motifs" << std::endl;
        std::cerr << "[" << name << "] counted " << count << " valid motifs in " << elapsed << "s" << std::endl;
    } else {
        std::cerr << "wrong type given, please provide one of these: 'AB','AF' or '-'" << std::endl;
    }
  }
  if (mode == 0) std::cerr << "total motifs counted: " << totalCount << std::endl;
}

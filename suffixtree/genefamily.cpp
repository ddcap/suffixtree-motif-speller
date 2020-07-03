
#include <string>
#include <vector>
#include <fstream>
#include "genefamily.h"

std::chrono::time_point<std::chrono::system_clock> prevTime;

void startChrono()
{
        prevTime = std::chrono::system_clock::now();
}

double stopChrono()
{
        std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - prevTime;
        return (elapsed.count());
}

void GeneFamily::readOrthologousFamily(const std::string& filename, std::pair<short, short> l, int maxDegeneration) {
    std::ifstream ifs(filename.c_str());
    readOrthologousFamily(ifs, l, maxDegeneration);
}

void GeneFamily::readOrthologousFamily(std::istream& ifs, std::pair<short, short> l, int maxDegeneration) {
  int totalCount = 0;
  while (ifs) {
    // READ DATA
    std::string T, newick, line, name;
    int N;
    getline(ifs, line);
    if(line.empty()) {continue;}
    name = line;
    getline(ifs, newick);
    getline(ifs, line);
    N = std::stoi(line);
    for (int i = 0; i < N; i++) {
        getline(ifs, line);
        getline(ifs, line);
        if (!T.empty())
            T.push_back('#');
        std::for_each(line.begin(), line.end(), [](char & c) { // convert all to upper case!
            c = ::toupper(c);
        });
        T.append(line);
        T.push_back('#');
        T.append(Motif::ReverseComplement(line));
    }
    T.push_back('$');
    // std::cerr << "T: " << T << std::endl;
    std::cerr << "[" << name << "] " << N << " gene families" << std::endl;
    // PROCESS DATA
    startChrono();
    BLSScore bls(newick);
    SuffixTree ST(T, true);
    int count = ST.printMotifsWithPositions(l, TWOFOLDSANDN, maxDegeneration, bls, std::cout);
    totalCount += count;
    double elapsed = stopChrono();
    std::cerr << "[" << name << "] counted " << count << " motifs in " << elapsed << "s" << std::endl;
  }
  std::cerr << "total motifs iterated: " << totalCount << std::endl;
}

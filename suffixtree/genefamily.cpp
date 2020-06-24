
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

void GeneFamily::readOrthologousFamily(const std::string& filename) {
    std::ifstream ifs(filename.c_str());
    readOrthologousFamily(ifs);
}

void GeneFamily::readOrthologousFamily(std::istream& ifs) {
  std::string T, newick, line, name;
  int N;
  std::pair<short, short> l(6, 13);
  int maxDegeneration = 3;
  while (ifs) {
    // READ DATA
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
        T.append(line);
        T.push_back('#');
        T.append(Motif::ReverseComplement(line));
    }
    T.push_back('$');
    // PROCESS DATA
    startChrono();
    BLSScore bls(newick);
    SuffixTree ST(T, true);
    ST.printMotifs(l, TWOFOLDSANDN, maxDegeneration, bls, std::cout);
    double elapsed = stopChrono();
    std::cerr << "time for iterating motifs: " << elapsed << std::endl;
  }
}

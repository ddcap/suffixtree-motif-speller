
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

void GeneFamily::readOrthologousFamily(const std::string& filename, bool typeIsAB, std::pair<short, short> l, int maxDegeneration) {
    std::ifstream ifs(filename.c_str());
    readOrthologousFamily(ifs, typeIsAB, l, maxDegeneration);
}

void GeneFamily::readOrthologousFamily(std::istream& ifs, bool typeIsAB, std::pair<short, short> l, int maxDegeneration) {
  int totalCount = 0;
  while (ifs) {
    std::vector<size_t> stringStartPositions;
    stringStartPositions.push_back(0);
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
            T.push_back(IupacMask::DELIMITER);
        std::for_each(line.begin(), line.end(), [](char & c) { // convert all to upper case!
            if(c == IupacMask::FILLER)
                c = IupacMask::DELIMITER;
            else
                c = ::toupper(c);
        });
        T.append(line);
        T.push_back(IupacMask::DELIMITER);
        stringStartPositions.push_back(T.size());
        T.append(Motif::ReverseComplement(line));
        stringStartPositions.push_back(T.size() + 1);
        // std::cout << T << std::endl;
    }
    T.push_back(IupacMask::DELIMITER);

    // std::cerr << "T: " << T << std::endl;
    std::cerr << "[" << name << "] " << N << " gene families" << std::endl;
    // PROCESS DATA
    startChrono();
    BLSScore bls(newick, N);
    // std::cout << T << std::flush;
    SuffixTree ST(T, true, stringStartPositions);

// TEST WRONG MOTIFS...
    // std::string testMotif = "SCGYCN";
    // occurence_bits occurence;
    // std::vector<std::pair<int, int>> testpos = ST.matchIupacPatternWithPositions(testMotif, bls, 3, occurence);
    // std::cerr << "testing " << testMotif << " with best occurence: " << +occurence  << " and BLS: " << bls.getBLSScore(occurence) << std::endl;
    // for (auto p: testpos) {
        // std::cerr << "[" << p.first << ", " << p.second << "]: " << ST.printPosPair(p, testMotif.length()) << std::endl;
    // }
    // int count = 0;


// FINAL CODE
    int count = ST.printMotifs(l, TWOFOLDSANDN, maxDegeneration, bls, std::cout, typeIsAB);
    size_t iteratorcount = ST.getMotifsIteratedCount();

    totalCount += count;
    double elapsed = stopChrono();
    std::cerr << "\33[2K\r[" << name << "] iterated over " << iteratorcount << " motifs" << std::endl;
    std::cerr << "[" << name << "] counted " << count << " valid motifs in " << elapsed << "s" << std::endl;

  }
  std::cerr << "total motifs counted: " << totalCount << std::endl;
}


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

void GeneFamily::readOrthologousFamily(std::ifstream& ifs) {
  std::string T, newick, line;
  int N;
  while (ifs) {
    // read the name
    getline(ifs, line);
    if(line.empty()) {continue;}
    // std::cerr << "processing "<< line;
    getline(ifs, newick);
    getline(ifs, line);
    N = std::stoi(line);
    // std::cerr << " with " << N << " families" << std::endl;
    for (int i = 0; i < N; i++) {
      getline(ifs, line);
      // std::cerr << line << " ";
      getline(ifs, line);
      if (!T.empty())
        T.push_back('#');
      T.append(line);
      T.push_back('#');
      T.append(Motif::ReverseComplement(line));
      // std::cerr << line << std::endl;
      // std::cerr << Motif::ReverseComplement(line) << std::endl;
      // std::cerr << i << "[" << (T.length() - (line.length()*2+1)) << ", " << (T.length() - (line.length()+1)) << "], RC[" << (T.length() - (line.length())) << ", " << T.length() << "]" << std::endl;
    }
  }
  T.push_back('$');

  // std::cerr << T << std::endl;
  // std::string teststr = "ACGTAGA";
  // std::cerr << teststr << " " << Motif::ReverseComplement(teststr) << std::endl;

// test NEWICK
  // std::cout << newick << std::endl;
  BLSScore bls(newick);
  // std::cerr << bls << std::endl;
  // std::bitset<N_BITS> test(0);
  // test.set(0);
  // std::cerr << "BLS of " << test << " is " << bls.getBLSScore(test) << std::endl;
  // test.set(1);
  // std::cerr << "BLS of " << test << " is " << bls.getBLSScore(test) << std::endl;
  // test.set(2);
  // std::cerr << "BLS of " << test << " is " << bls.getBLSScore(test) << std::endl;
  // test.set(3);
  // std::cerr << "BLS of " << test << " is " << bls.getBLSScore(test) << std::endl;
  // test.reset(); test.set(3); test.set(1);
  // std::cerr << "BLS of " << test << " is " << bls.getBLSScore(test) << std::endl;
  // test.set(2);
  // std::cerr << "BLS of " << test << " is " << bls.getBLSScore(test) << std::endl;
  // int count = 999;
  // std::vector<std::bitset<N_BITS>> testset;
  // for(int i =0; i < count; i++) {
  //     testset.push_back(bls.random_bitset(4));
  // }
  // startChrono();
  // for(int j =0; j < 10; j++) {
  //   for(int i =0; i < count; i++) {
  //       bls.getBLSScore(testset[i]);
  //   }
  // }
  // double elapsed = stopChrono();
  // std::cout << "Time for " << count << " bls scores: " << elapsed << std::endl;


// create ST
  // std::cerr << "Building suffix tree..." << std::endl;
    SuffixTree ST(T, true);
  // std::string word = "AAAATCTTGTTT";
  // ST.matchPattern(word, bls);

    // std::string iupacword = "AAAAKMTKGTTT";
    // std::bitset<N_BITS> occurence(0);
    // std::unordered_set<STPosition, STPositionHash> positions;
    // ST.matchIupacPattern(iupacword, positions, occurence);
    //
    // if(!positions.empty()) {
    // std::cerr << iupacword << " found with occurence " << occurence << " at " << positions.size() << " locations" << std::endl;
    // for(auto p : positions) {
    //     std::cerr << "pos: " << p.getPositionInText() << ": " <<  T.substr(p.getPositionInText() - p.getDepth(), iupacword.length()) << std::endl;
    // }
    // }
    // std::string iupacword2 = "AAACMAKMTTTT";
    // ST.matchIupacPattern(iupacword2, positions, occurence);
    //
    // if(!positions.empty()) {
    //   std::cerr << iupacword2 << " found with occurence " << occurence << " at " << positions.size() << " locations" << std::endl;
    //   for(auto p : positions) {
    //       std::cerr << "pos: " << p.getPositionInText() << ": " <<  T.substr(p.getPositionInText() - p.getDepth(), iupacword.length()) << std::endl;
    //   }
    // }

  // std::cerr << ST << std::endl;
  std::pair<short, short> l(6, 13);
  int maxDegeneration = 3;
  float minBlsScore = 0.15; // min threshold
  std::ofstream myfile;
  myfile.open ("/tmp/motifs.txt");
  ST.printMotifs(l, TWOFOLDSANDN, maxDegeneration, bls, minBlsScore, myfile);
  myfile.close();

}


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
  std::string T, newick, line, name;
  int N;
  while (ifs) {
    // read the name
    getline(ifs, line);
    if(line.empty()) {continue;}
    name = line;
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
  // for(int j =0; j < 10; j++) {
  //   for(int i =0; i < count; i++) {
  //       bls.getBLSScore(testset[i]);
  //   }
  // }




// create ST
  // std::cerr << "Building suffix tree..." << std::endl;
    startChrono();
    SuffixTree ST(T, true);
  // std::string word = "AAAATCTTGTTT";
  // ST.matchPattern(word, bls);

    // std::string iupacword = "TNAGCN";
    // occurence_bits occurence(0);
    // std::vector<STPosition> positions = ST.matchIupacPattern(iupacword, 3, occurence);
    //
    // if(!positions.empty()) {
    //     std::cerr << iupacword << " found with occurence " << occurence << " [" << bls.getBLSScore(occurence) << "] at " << positions.size() << " locations" << std::endl;
    //     for(auto p : positions) {
    //         std::cerr << "pos: " << p.getPositionInText() << ": " <<  T.substr(p.getPositionInText() - p.getDepth(), iupacword.length()) << std::endl;
    //     }
    // }
    // std::string iupacword2 = "AAACMAKNTTTT";
    // ST.matchIupacPattern(iupacword2, positions, occurence);
    //
    // if(!positions.empty()) {
    //   std::cerr << iupacword2 << " found with occurence " << occurence << " [" << bls.getBLSScore(occurence) << "] at " << positions.size() << " locations" << std::endl;
    //   for(auto p : positions) {
    //       std::cerr << "pos: " << p.getPositionInText() << ": " <<  T.substr(p.getPositionInText() - p.getDepth(), iupacword.length()) << std::endl;
    //   }
    // }
    // std::string iupacword3 = "AAAANMTKGTTT";
    // ST.matchIupacPattern(iupacword3, positions, occurence);
    //
    // if(!positions.empty()) {
    //   std::cerr << iupacword3 << " found with occurence " << occurence << " [" << bls.getBLSScore(occurence) << "] at " << positions.size() << " locations" << std::endl;
    //   for(auto p : positions) {
    //       std::cerr << "pos: " << p.getPositionInText() << ": " <<  T.substr(p.getPositionInText() - p.getDepth(), iupacword.length()) << std::endl;
    //   }
    // }

  // std::cerr << ST << std::endl;
  std::pair<short, short> l(6, 13);
  int maxDegeneration = 3;
  // std::ofstream myfile;
  // myfile.open ("/tmp/motifs.txt");
  ST.printMotifs(l, TWOFOLDSANDN, maxDegeneration, bls, std::cout);
  // myfile.close();
  double elapsed = stopChrono();
  std::cerr << "time for iterating motifs: " << elapsed << std::endl;

}

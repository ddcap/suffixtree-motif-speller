#ifndef GENEFAMILY_H
#define GENEFAMILY_H


#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include "suffixtree.h"


#define MAX_VALID_CHARS 5
class GeneFamily {
private:

    static const std::unordered_set<char> validCharacters;

public:
  static void readOrthologousFamily(const std::string& filename, std::vector<float> blsThresholds_, Alphabet alphabet, bool typeIsAB, std::pair<short, short> l, int maxDegeneration);

  static void readOrthologousFamily(std::istream& ifs, std::vector<float> blsThresholds_, Alphabet alphabet, bool typeIsAB, std::pair<short, short> l, int maxDegeneration) ;
};

#endif

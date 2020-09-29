#ifndef GENEFAMILY_H
#define GENEFAMILY_H


#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include "suffixtree.h"


class GeneFamily {
private:

public:
  static void readOrthologousFamily(const std::string& filename, std::vector<float> blsThresholds_, Alphabet alphabet, bool typeIsAB, std::pair<short, short> l, int maxDegeneration);

  static void readOrthologousFamily(std::istream& ifs, std::vector<float> blsThresholds_, Alphabet alphabet, bool typeIsAB, std::pair<short, short> l, int maxDegeneration) ;
};

#endif

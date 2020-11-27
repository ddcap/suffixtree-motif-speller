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

    static size_t getIndexOfVector(const std::vector<std::string> &v, const std::string &val);
public:
    static void readOrthologousFamily(const int mode, const std::string& filename, const std::vector<float> blsThresholds_, const Alphabet alphabet, const int type, const std::pair<short, short> l, const int maxDegeneration);

    static void readOrthologousFamily(const int mode, std::istream& ifs, const std::vector<float> blsThresholds_, const Alphabet alphabet, const int type, const std::pair<short, short> l, const int maxDegeneration);

    static void readGenes(std::istream& ifs, const int maxDegeneration, const short maxLen);
};

#endif

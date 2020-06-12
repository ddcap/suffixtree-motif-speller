#ifndef GENEFAMILY_H
#define GENEFAMILY_H


#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include "suffixtree.h"


class GeneFamily {
private:
    static const std::vector<unsigned char> complement;

    static std::string RC(std::string read);

public:
  static void readOrthologousFamily(const std::string& filename);

  static void readOrthologousFamily(std::ifstream& ifs) ;
};

#endif

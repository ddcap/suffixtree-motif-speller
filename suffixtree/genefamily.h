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
  static void readOrthologousFamily(const std::string& filename);

  static void readOrthologousFamily(std::istream& ifs) ;
};

#endif

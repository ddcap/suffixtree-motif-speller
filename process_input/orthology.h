#ifndef ORTHOLOGY_H
#define ORTHOLOGY_H


#include <iostream>
#include <fstream>
#include <list>
#include <stack>
#include "genes.h"
#include "newick.h"

class Orthology {
private:
  const char delim = ' ';
  void formatGeneList(std::ostream& o, Genes *genemap, Newick *newick, std::string cluster, std::string genelist);
  void readOrthology(std::ostream& o, Genes *genes, Newick *newick, std::string orthologyFile);
public:
  Orthology(Genes *genes, Newick *newick, std::string orthologyFile) {
    readOrthology(std::cout, genes, newick, orthologyFile);
  }
};


#endif

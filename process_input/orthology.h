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
  void readOrthology(Genes *genemap, Newick *newick, std::string orthologyFile, std::string outputfolder);
public:
  Orthology(Genes *genes, Newick *newick, std::string orthologyFile, std::string outputfolder) {
    readOrthology(genes, newick, orthologyFile, outputfolder);
  }
};


#endif

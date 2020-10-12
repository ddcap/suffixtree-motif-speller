#ifndef GENES_H
#define GENES_H


#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <stack>
#include <vector>

struct Gene {
  std::string species;
  std::string sequence;
  friend std::ostream& operator<< (std::ostream& o, const Gene& b) {
    o << b.species << "\t" << b.sequence;
    return o;
  }
  bool operator==(const Gene &other) const {
    return (species == other.species && sequence == other.sequence);
  }
};
namespace std {
  template <>
  struct hash<Gene>
  {
    std::size_t operator()(const Gene& k) const
    {
      using std::hash;
      using std::string;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:

      return ((hash<string>()(k.species)
               ^ (hash<string>()(k.sequence) << 1)) >> 1);
    }
  };
}


class Genes {
private:
  static const std::vector<std::string> characterToMask;
  static const std::unordered_set<char> validCharacters;
  std::unordered_map<std::string, Gene> genemap;

  static void fixSequence(std::string &seq);
  static char getRandomCharFromIupac(char c);
  void readFastas(std::string directory);
  void readFasta(std::string species, std::string fasta);
public:
  Genes(std::string directory) {
    readFastas(directory);
  }
  Gene *getGene(std::string geneid);
};


#endif

#ifndef MOTIFMAP_H
#define MOTIFMAMP_H

#include "motif.h"

#define IUPAC_FULL_COUNT 15

class MotifMap { // first version
private:
  MotifMap *children;
  blscounttype *v;
  void createBlsVectorFromByte(const int& val,const char &blsvectorsize);
  void addByteToBlsVector(const int& val);
  void initChildren() { children = new MotifMap[IUPAC_FULL_COUNT]; }

public:
  MotifMap(): children(NULL), v(NULL) {}
  MotifMap(const char &blsvectorsize): children(NULL), v(NULL) {}
  void addMotifToMap(const std::string &motif, const size_t pos, const int &val,const char &blsvectorsize);
  void recPrintAndDelete(const std::string currentmotif, long &unique_count, std::ostream &out, const short &maxlen, const char &blsvectorsize);
};


class SparseMotifMapNode {
private:
  SparseMotifMapNode *data; // Use a single pointer to reduce memory usage, since we get millions of these very small objects, so 8byte vs 16/24 bytes matters a lot
  // DATA ordered like this:
  // first few bytes are iupacchar to index map
  // next few bytes are the blsvector if there is one (can be 0) -> index.first == index.second
  // next bytes are actual pointers to children
  // ------------------------------------------------------------------------
  // |iupac to index map(children)|bls vector of this node|children pointers|
  // ------------------------------------------------------------------------
  void addByteToBlsVector(const int& val, const std::pair<int, int> &startIndexes);
  void addChild(const int &iupac_value, const std::pair<int, int> &startIndexes, const int &pos, const std::pair<short, short> &range, const char &blsvectorsize);
  void init(const std::pair<int, int> &startIndexes);

public:
  ~SparseMotifMapNode() { free(data); };
  SparseMotifMapNode(const std::pair<int, int> &startIndexes);
  SparseMotifMapNode(const int &startIndex);
  void addMotifToMap(const std::string &motif, const size_t pos, const int &val, const std::pair<int, int> &startIndexes, const std::pair<short, short> &range, const char &blsvectorsize);
  void recPrintAndDelete(const std::string currentmotif, long &unique_count, std::ostream &out, const std::pair<int, int> &startIndexes, const std::pair<short, short> &range, const char &blsvectorsize);
};
// create another class that also has its own bls vector! for if range allows multiple lengths!
class MotifMapLeafs {
private:
  char *data;
  // DATA ordered like this:
  // first char 0 -> no bls vector for this node, 1 this node also has a bls vector -> is first one in list
  // first few bytes are iupacchar to index map
  // next bytes are actual blsvector ordered according to previous map!
  // --------------------------------------------------------------------
  // |has bls vector bool|iupac to index map|bls vectors in order of map|
  // --------------------------------------------------------------------
  void createNewBlsVector(const int& iupac_value, const int& val, const char &blsvectorsize);
public:
  void init();
  MotifMapLeafs() : data(NULL) { init(); };
  void addToBlsVector(const int& iupac_value, const int& val, const char &blsvectorsize);
  void printMotifsAndDeleteData(const std::string currentmotif, long &unique_count, std::ostream &out, const std::pair<short, short> &range, const char &blsvectorsize);
};

class SparseMotifMap { // second version, less memory usage
private:
  SparseMotifMapNode *root;
  const char blsvectorsize;
  const std::pair<int, int> startIndexes; // first index is start index of blsvector, second index is start index of actual nodes
  const std::pair<short, short> range;
public:
  SparseMotifMap(const char &blsvectorsize, const std::pair<short, short> &range);
  void addMotifToMap(const std::string &motif, const int &val);
  void recPrintAndDelete(long &unique_count, std::ostream &out);
};
#endif

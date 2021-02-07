#ifndef MOTIFMAP_H
#define MOTIFMAMP_H

#include "motif.h"

#define IUPAC_FULL_COUNT 15

class MotifMap {
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
  SparseMotifMapNode *data;
  // DATA ordered like this:
  // first few bytes are iupacchar to index map
  // next few bytes are the blsvector if there is one (can be 0) -> index.first == index.second
  // next bytes are actual pointers to children
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
// create another class that also has its own bls vector! for if range allows allows multiple lengths!
class MotifMapLeafs {
private:
  char *data;
  // DATA ordered like this:
  // first char 0 -> no bls vector for this node, 1 this node also has a bls vector -> is first one in list
  // first few bytes are iupacchar to index map
  // next bytes are actual blsvector ordered according to previous map!
  void createNewBlsVector(const int& iupac_value, const int& val, const char &blsvectorsize);
  void freeData() { free(data); }
public:
  void init();
  MotifMapLeafs() : data(NULL) { init(); };
  ~MotifMapLeafs() { free(data); };
  void addToBlsVector(const int& iupac_value, const int& val, const char &blsvectorsize);
  void printMotifsAndDeleteData(const std::string currentmotif, long &unique_count, std::ostream &out, const std::pair<short, short> &range, const char &blsvectorsize);
};

class SparseMotifMap {
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
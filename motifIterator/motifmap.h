#ifndef MOTIFMAP_H
#define MOTIFMAMP_H

#include "motif.h"

#define IUPAC_FULL_COUNT 15
class MotifMap {
private:
  MotifMap *children;
  blscounttype *v;
  void createBlsVectorFromByte(const int& val, const int &listSize);
  void addByteToBlsVector(const int& val);
  void initChildren() { children = new MotifMap[IUPAC_FULL_COUNT]; }

public:
  MotifMap(): children(NULL), v(NULL) {}
  void addMotifToMap(const std::string &motif, const size_t pos, const int &val, const int &listSize);
  void recPrintAndDelete(const std::string currentmotif, long &unique_count, std::ostream &out, const short &maxlen, const char &blsvectorsize);
};

class SparseMotifMap {
private:
  SparseMotifMap *children;
  // char *childrenChars; // first element here is the size! so we have 1/4 less memory per node!
  // blscounttype *v;
  void createBlsVectorFromByte(const int& val, const char &blsvectorsize);
  void addByteToBlsVector(const int& val, const char &blsvectorsize);
  void addChild(int iupac_value, const char &blsvectorsize);
  int getFirstIndexOfChild(const char &blsvectorsize);
  SparseMotifMap *getChildren(const char &blsvectorsize);
  char *getBlsVector(const char &blsvectorsize);
  char *getIupacMapping();
public:
  SparseMotifMap(const char &blsvectorsize);
  void addMotifToMap(const std::string &motif, const size_t pos, const int &val, const char &blsvectorsize);
  void recPrintAndDelete(const std::string currentmotif, long &unique_count, std::ostream &out, const short &maxlen, const char &blsvectorsize);
};

#endif

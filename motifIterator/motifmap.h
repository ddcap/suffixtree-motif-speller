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
  char *childrenChars; //
  unsigned char childrencount;
  blscounttype *v;
  void createBlsVectorFromByte(const int& val, const int &listSize);
  void addByteToBlsVector(const int& val);
  void initChildren(int iupac_value);
  void addChild(int iupac_value);

public:
  SparseMotifMap(): children(NULL), childrenChars(NULL), childrencount(0), v(NULL) { }
  void addMotifToMap(const std::string &motif, const size_t pos, const int &val, const int &listSize);
  void recPrintAndDelete(const std::string currentmotif, long &unique_count, std::ostream &out, const short &maxlen, const char &blsvectorsize);
};

#endif

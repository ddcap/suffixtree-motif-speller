#include "motifmap.h"

// MOTIFMAP
void MotifMap::addMotifToMap(const std::string &motif, const size_t pos, const int &val, const int &listSize) {
  if(pos == motif.length()) {
    if(v == NULL) {
      createBlsVectorFromByte(val, listSize);
    } else {
      addByteToBlsVector(val);
    }
  } else {
    if(children == NULL) {
      initChildren();
    }
    children[IupacMask::characterToMask[motif[pos]].getMask() - 1].addMotifToMap(motif, pos + 1, val, listSize);
  }
}


void MotifMap::recPrintAndDelete(const std::string currentmotif, long &unique_count, std::ostream &out, const short &maxlen, const char &blsvectorsize) {
    if(v != NULL) {
        Motif::writeGroupIDAndMotifInBinary(currentmotif, maxlen, out);
        out.write(&blsvectorsize, 1); // assume
        for(int i = 0; i < blsvectorsize; i++) {
            std::cout.write((char*)&v[i], sizeof(blscounttype));
        }
        delete v;
        unique_count++;
    }
    if(children != NULL) {
      for (int i = 0; i < IUPAC_FULL_COUNT; i++) {
        children[i].recPrintAndDelete(currentmotif + IupacMask::representation[i + 1], unique_count, out, maxlen, blsvectorsize);
      }
      delete children;
    }
}

void MotifMap::createBlsVectorFromByte(const int& val, const int &listSize) {
  v = new blscounttype[listSize];
  int i = 0;
  for (i = 0; i < val; i++) {
      v[i] = 1;
  }
  for (; i < listSize; i++) {
      v[i] = 0; // make sure its 0
  }
}

void MotifMap::addByteToBlsVector(const int& val) {
    for (int i = 0; i < val; i++) {
        v[i] += 1;
    }
}

// SPARSEMOTIFMAP
void SparseMotifMap::addMotifToMap(const std::string &motif, const size_t pos, const int &val, const int &listSize) {
  if(pos == motif.length()) {
    if(v == NULL) {
      createBlsVectorFromByte(val, listSize);
    } else {
      addByteToBlsVector(val);
    }
  } else {
    int iupac_value = IupacMask::characterToMask[motif[pos]].getMask() - 1;
    if(children == NULL) {
      initChildren(iupac_value);
    } else if(childrenChars[iupac_value] == -1) {
      addChild(iupac_value);
    }
    children[(short)childrenChars[iupac_value]].addMotifToMap(motif, pos + 1, val, listSize);
  }
}


void SparseMotifMap::recPrintAndDelete(const std::string currentmotif, long &unique_count, std::ostream &out, const short &maxlen, const char &blsvectorsize) {
    if(v != NULL) {
        Motif::writeGroupIDAndMotifInBinary(currentmotif, maxlen, out);
        out.write(&blsvectorsize, 1); // assume
        for(int i = 0; i < blsvectorsize; i++) {
            std::cout.write((char*)&v[i], sizeof(blscounttype));
        }
        delete v;
        unique_count++;
    }
    if(children != NULL) {
      for (int i = 0; i < IUPAC_FULL_COUNT; i++) {
          if(childrenChars[i] != -1) {
              children[childrenChars[i]].recPrintAndDelete(currentmotif + IupacMask::representation[i + 1], unique_count, out, maxlen, blsvectorsize);
          }
      }
      delete children;
      delete childrenChars;
    }
}

void SparseMotifMap::initChildren(int iupac_value) {
    children = new SparseMotifMap[1];
    childrenChars = new char[IUPAC_FULL_COUNT]{ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    childrenChars[iupac_value] = childrencount;
    childrencount++;
}

void SparseMotifMap::addChild(int iupac_value) {
    childrenChars[iupac_value] = childrencount;
    childrencount++;
    children = (SparseMotifMap *)realloc(children, childrencount * sizeof (SparseMotifMap));
    new ((void*)&children[childrencount - 1]) SparseMotifMap(); // call constructor
}


void SparseMotifMap::createBlsVectorFromByte(const int& val, const int &listSize) {
    v = new blscounttype[listSize];
    int i = 0;
    for (i = 0; i < val; i++) {
        v[i] = 1;
    }
    for (; i < listSize; i++) {
        v[i] = 0; // make sure its 0
    }
}

void SparseMotifMap::addByteToBlsVector(const int& val) {
    for (int i = 0; i < val; i++) {
        v[i] += 1;
    }
}

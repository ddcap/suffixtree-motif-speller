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
SparseMotifMap::SparseMotifMap(const char &blsvectorsize) { //, childrenChars(NULL), v(NULL) { }
    int mappingsize_in_SpareMotifMaps = 1 + ((sizeof(char) * (IUPAC_FULL_COUNT + 1)) - 1) / sizeof(SparseMotifMap);
    int blsvector_in_SpareMotifMaps = 1 + ((sizeof(blscounttype) * blsvectorsize) - 1) /  sizeof(SparseMotifMap);
    children = (SparseMotifMap *)malloc((mappingsize_in_SpareMotifMaps + blsvector_in_SpareMotifMaps) * sizeof(SparseMotifMap));
    memset(children, 0, (mappingsize_in_SpareMotifMaps + blsvector_in_SpareMotifMaps) * sizeof(SparseMotifMap));
    memset(&((char *)&children[0])[1], -1, IUPAC_FULL_COUNT);
    // char *iupac_mapping = getIupacMapping();
    // char *v = getBlsVector(blsvectorsize);
    // for (int i = 0 ; i < blsvectorsize; i++ ) {
    //     std::cerr << "v[" << i << "]: " << +v[i] << std::endl;
    // }
    // for (int i = 0 ; i < IUPAC_FULL_COUNT + 1; i++ ) {
    //     std::cerr << "iupac_mapping[" << i << "]: " << +iupac_mapping[i] << std::endl;
    // }
}

SparseMotifMap *SparseMotifMap::getChildren(const char &blsvectorsize) {
    return (SparseMotifMap *)&children[1 + ((sizeof(char) * (IUPAC_FULL_COUNT + 1)) - 1) / sizeof(SparseMotifMap) + 1 + ((sizeof(blscounttype) * blsvectorsize) - 1) /  sizeof(SparseMotifMap)];
}
char *SparseMotifMap::getBlsVector(const char &blsvectorsize) {
    return (char *)&children[1 + ((sizeof(char) * (IUPAC_FULL_COUNT + 1)) - 1) / sizeof(SparseMotifMap)];
}
char *SparseMotifMap::getIupacMapping() {
    return (char *)&children[0];
}

void SparseMotifMap::addMotifToMap(const std::string &motif, const size_t pos, const int &val, const char &blsvectorsize) {
    if(pos == motif.length()) {
        if(getBlsVector(blsvectorsize)[0] == -1) {
            createBlsVectorFromByte(val, blsvectorsize);
        } else {
            addByteToBlsVector(val, blsvectorsize);
        }
    } else {
        int iupac_value = IupacMask::characterToMask[motif[pos]].getMask();
        // std::cerr << "iupac_val " << iupac_value << std::endl;
        if(getIupacMapping()[iupac_value] == -1) {
            addChild(iupac_value, blsvectorsize);
            // std::cerr << "new iupac_mapping; " << (short)getIupacMapping()[iupac_value] << std::endl;
        }
        // need to do getIupacMapping() again cause the address might have changed!!!
        int startIdx = 1 + ((sizeof(char) * (IUPAC_FULL_COUNT + 1)) - 1) / sizeof(SparseMotifMap) + 1 + ((sizeof(blscounttype) * blsvectorsize) - 1) /  sizeof(SparseMotifMap);
        children[startIdx + (short)getIupacMapping()[iupac_value]].addMotifToMap(motif, pos + 1, val, blsvectorsize);
    }
}


void SparseMotifMap::recPrintAndDelete(const std::string currentmotif, long &unique_count, std::ostream &out, const short &maxlen, const char &blsvectorsize) {
    int startIdx = 1 + ((sizeof(char) * (IUPAC_FULL_COUNT + 1)) - 1) / sizeof(SparseMotifMap) + 1 + ((sizeof(blscounttype) * blsvectorsize) - 1) /  sizeof(SparseMotifMap);
    char *iupac_mapping = getIupacMapping();
    char *v = getBlsVector(blsvectorsize);
    if(getBlsVector(blsvectorsize)[0] > 0) {
        Motif::writeGroupIDAndMotifInBinary(currentmotif, maxlen, out);
        out.write(&blsvectorsize, 1); // assume
        for(int i = 0; i < blsvectorsize; i++) {
            std::cout.write((char*)&v[i], sizeof(blscounttype));
        }
        unique_count++;
    }
    if(children != NULL) {
      for (int i = 0; i < IUPAC_FULL_COUNT; i++) {
          if(iupac_mapping[i + 1] != -1) {
              children[startIdx + iupac_mapping[i + 1]].recPrintAndDelete(currentmotif + IupacMask::representation[i + 1], unique_count, out, maxlen, blsvectorsize);
          }
      }
      delete children;
    }
}

void SparseMotifMap::addChild(int iupac_value, const char &blsvectorsize) {
    char *iupac_mapping = getIupacMapping();
    int startIdx = 1 + ((sizeof(char) * (IUPAC_FULL_COUNT + 1)) - 1) / sizeof(SparseMotifMap) + 1 + ((sizeof(blscounttype) * blsvectorsize) - 1) /  sizeof(SparseMotifMap);
    iupac_mapping[iupac_value] = iupac_mapping[0];
    iupac_mapping[0]++;
    children = (SparseMotifMap *)realloc(children, (startIdx + getIupacMapping()[0]) * sizeof (SparseMotifMap));
    new ((void*)&children[startIdx + getIupacMapping()[0] - 1]) SparseMotifMap(blsvectorsize); // call constructor
}


void SparseMotifMap::createBlsVectorFromByte(const int& val, const char &blsvectorsize) {
    char *v = getBlsVector(blsvectorsize);
    int i = 0;
    for (i = 0; i < val; i++) {
        v[i] = 1;
    }
    for (; i < blsvectorsize; i++) {
        v[i] = 0; // make sure its 0
    }
    // for (int i = 0 ; i < blsvectorsize; i++ ) {
    //     std::cerr << "v[" << i << "]: " << +v[i] << std::endl;
    // }
}

void SparseMotifMap::addByteToBlsVector(const int& val, const char &blsvectorsize) {
    char *v = getBlsVector(blsvectorsize);
    for (int i = 0; i < val; i++) {
        v[i] += 1;
    }
    // for (int i = 0 ; i < blsvectorsize; i++ ) {
    //     std::cerr << "v[" << i << "]: " << +v[i] << std::endl;
    // }
}

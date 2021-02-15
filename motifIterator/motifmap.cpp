#include "motifmap.h"
#include "malloc.h"

// MOTIFMAP
void MotifMap::addMotifToMap(const std::string &motif, const size_t pos, const int &val,const char &blsvectorsize) {
  if(pos == motif.length()) {
    if(v == NULL) {
      createBlsVectorFromByte(val, blsvectorsize);
    } else {
      addByteToBlsVector(val);
    }
  } else {
    if(children == NULL) {
      initChildren();
    }
    children[IupacMask::characterToMask[motif[pos]].getMask() - 1].addMotifToMap(motif, pos + 1, val, blsvectorsize);
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

void MotifMap::createBlsVectorFromByte(const int& val, const char &blsvectorsize) {
  v = new blscounttype[blsvectorsize];
  int i = 0;
  for (i = 0; i < val; i++) {
      v[i] = 1;
  }
  for (; i < blsvectorsize; i++) {
      v[i] = 0; // make sure its 0
  }
}

void MotifMap::addByteToBlsVector(const int& val) {
    for (int i = 0; i < val; i++) {
        v[i] += 1;
    }
}

// SPARSEMOTIFMAPNODE
SparseMotifMapNode::SparseMotifMapNode(const std::pair<int, int> &startIndexes) {
    init(startIndexes);
}
SparseMotifMapNode::SparseMotifMapNode(const int &startIndex) {
    init(std::pair<int, int>(startIndex, startIndex));
}

void SparseMotifMapNode::init(const std::pair<int, int> &startIndexes) {
    data = (SparseMotifMapNode *)malloc((startIndexes.second) * sizeof(SparseMotifMap));
    memset(data, 0, (startIndexes.second) * sizeof(SparseMotifMap));
    memset(&((char *)&data[0])[1], -1, IUPAC_FULL_COUNT); // set map to -1 which means this iupac letter isnt set yet
}
void SparseMotifMapNode::addByteToBlsVector(const int& val, const std::pair<int, int> &startIndexes) {
    blscounttype *v = (blscounttype *)&data[startIndexes.first];
    for (int i = 0; i < val; i++) {
        v[i] += 1;
    }
}

void SparseMotifMapNode::addChild(const int &iupac_value, const std::pair<int, int> &startIndexes, const int &pos, const std::pair<short, short> &range, const char &blsvectorsize) {
    char *iupac_mapping = (char *)&data[0];
    iupac_mapping[iupac_value] = iupac_mapping[0];
    iupac_mapping[0]++;
    data = (SparseMotifMapNode *)realloc(data, (startIndexes.second + iupac_mapping[0]) * sizeof (SparseMotifMapNode));
}

void SparseMotifMapNode::addMotifToMap(const std::string &motif, const size_t pos, const int &val, const std::pair<int, int> &startIndexes, const std::pair<short, short> &range, const char &blsvectorsize) {
    char *iupac_mapping = (char *)&data[0];
    int iupac_value = IupacMask::characterToMask[motif[pos]].getMask();
    // std::cerr << motif << "-> " << pos << " vs " << range.first  << " - " << range.second  << " position in idx map: " << +iupac_mapping[iupac_value]<< std::endl;
    if(iupac_mapping[iupac_value] == -1) {
        addChild(iupac_value, startIndexes, pos, range, blsvectorsize);
        iupac_mapping = (char *)&data[0];  // realloced so need a new address!
        if(pos < range.first - 1 && pos < range.second - 3) { // TODO fix this logic for ranges with more then 1 motif length!!!
            if (pos < range.first - 1) { // under the range of valid motifs
                // std::cerr << "creating node for " << motif.substr(0, pos + 1) << std::endl;
                data[startIndexes.second + iupac_mapping[iupac_value]].init(std::pair<int, int>(startIndexes.first, startIndexes.first));
            } else { // also needs a bls vector cause its in range of valid motifs
                // std::cerr << "creating node with vector for " << motif.substr(0, pos + 1) << std::endl;
                data[startIndexes.second + iupac_mapping[iupac_value]].init(startIndexes);
            }
        } else {  // end the range of valid motifs
            // std::cerr << "creating leafs node for " << motif.substr(0, pos + 1) << std::endl;
            ((MotifMapLeafs *)&data[startIndexes.second + iupac_mapping[iupac_value]])->init();
        }
    }
    if (pos < range.second - 3){ // in range of full nodes
        data[startIndexes.second + (short)iupac_mapping[iupac_value]].addMotifToMap(motif, pos + 1, val, startIndexes, range, blsvectorsize);
    } else {  // end the range of valid motifs
        // std::cerr <<"leafs node address for " << motif.substr(0, pos + 1) << ": " << &data[startIndexes.second + (short)iupac_mapping[iupac_value]] << std::endl;
        ((MotifMapLeafs *)&data[startIndexes.second + (short)iupac_mapping[iupac_value]])->addToBlsVector(IupacMask::characterToMask[motif[pos + 1]].getMask(), val, blsvectorsize);
    }
}


void SparseMotifMapNode::recPrintAndDelete(const std::string currentmotif, long &unique_count, std::ostream &out, const std::pair<int, int> &startIndexes, const std::pair<short, short> &range, const char &blsvectorsize) {
    char *iupac_mapping = (char *)&data[0];
    blscounttype *v = (blscounttype *)&data[startIndexes.first];
    size_t pos = currentmotif.length();
    if( pos < range.second - 3){
        if (pos >= range.first - 1  && pos < range.second - 3){ // bls vector inside this node
        // if(v[0] > 0) { // should always be the case though!
            Motif::writeGroupIDAndMotifInBinary(currentmotif, range.second, out);
            out.write(&blsvectorsize, 1); // assume
            for(int i = 0; i < blsvectorsize; i++) {
                std::cout.write((char*)&v[i], sizeof(blscounttype));
            }
            unique_count++;
        }
        for (int i = 0; i < IUPAC_FULL_COUNT; i++) {
            if(iupac_mapping[i + 1] != -1) {
                data[startIndexes.second + iupac_mapping[i + 1]].recPrintAndDelete(currentmotif + IupacMask::representation[i + 1], unique_count, out, startIndexes, range, blsvectorsize);
            }
        }
        free(data);
    } else {  // bls vector in motifmapleafs
        for (int i = 0; i < IUPAC_FULL_COUNT; i++) {
            if(iupac_mapping[i + 1] != -1) {
                ((MotifMapLeafs *)&data[startIndexes.second + iupac_mapping[i + 1]])->printMotifsAndDeleteData(currentmotif + IupacMask::representation[i + 1], unique_count, out, range, blsvectorsize);
            }
        }
        free(data); // WHY THIS NO GOOD????
    }
    if(pos == 2)  {
    //    std::cerr << "freeing memory of node " <<  currentmotif << std::endl;
        malloc_trim(0); // this gives memory back to OS!
    }
}

// MOTIFMAPLEAFS
void MotifMapLeafs::init() {
    data = (char *)malloc(1 + IUPAC_FULL_COUNT);
    data[0] = 0;
    memset(&((char *)&data[0])[1], -1, IUPAC_FULL_COUNT); // set map to -1 which means this iupac letter isnt set yet
}
void MotifMapLeafs::createNewBlsVector(const int& iupac_value, const int& val, const char &blsvectorsize) {
    this->data[iupac_value] = data[0];
    data[0]++;
    data = (char *)realloc(data, 1 + IUPAC_FULL_COUNT + blsvectorsize*sizeof(blscounttype)*data[0]);
    memset(&data[1 + IUPAC_FULL_COUNT + blsvectorsize*sizeof(blscounttype)*(data[0] - 1)], 0, blsvectorsize*sizeof(blscounttype)); // init bls vector to 0!
    // std::cerr << iupac_value << " new vector\n";
}

void MotifMapLeafs::addToBlsVector(const int& iupac_value, const int& val, const char &blsvectorsize) {
    char *iupac_mapping = (char *)&data[0];
    if(iupac_mapping[iupac_value] == -1) {
        createNewBlsVector(iupac_value, val, blsvectorsize);
        iupac_mapping = (char *)&data[0]; // realloced so need a new address!
    }
    // std::cerr << "current iupac[" << (void *)&data[0] << "] to idx map ";
    // for (int i = 0; i < IUPAC_FULL_COUNT; i++) {
    //     std::cerr << +iupac_mapping[i + 1] << " ";
    // }
    // std::cerr << "\n";
    blscounttype *v = (blscounttype *)&data[1 + IUPAC_FULL_COUNT + blsvectorsize*sizeof(blscounttype)*data[iupac_value]];
    // std::cerr << iupac_value << " added vector " << val << ": ";
    for (int i = 0; i < val; i++) {
        v[i] += 1;
        // std::cerr << +v[i] << " ";
    }
    // std::cerr << "\n";
}

void MotifMapLeafs::printMotifsAndDeleteData(const std::string currentmotif, long &unique_count, std::ostream &out, const std::pair<short, short> &range, const char &blsvectorsize) {
    char *iupac_mapping = (char *)&data[0];
    for (int i = 0; i < IUPAC_FULL_COUNT; i++) {
        if(iupac_mapping[i + 1] != -1) {
            blscounttype *v = (blscounttype *)&data[1 + IUPAC_FULL_COUNT + blsvectorsize*sizeof(blscounttype)*data[i + 1]];
            Motif::writeGroupIDAndMotifInBinary(currentmotif + IupacMask::representation[i + 1], range.second, out);
            out.write(&blsvectorsize, 1); // assume less than 256
            for(int i = 0; i < blsvectorsize; i++) {
                std::cout.write((char*)&v[i], sizeof(blscounttype));
            }
            unique_count++;
     //       if(unique_count % 1000000 == 0) {
     //           std::cerr << "counted " << unique_count << std::endl;
     //       }
        }
    }
    free(data);
}

// SPARSEMOTIFMAP
SparseMotifMap::SparseMotifMap(const char &blsvectorsize, const std::pair<short, short> &range) : root(NULL), blsvectorsize(blsvectorsize),
startIndexes(std::pair<int,int>(
    1 + ((sizeof(char) * (IUPAC_FULL_COUNT + 1)) - 1) / sizeof(SparseMotifMapNode), // index of bls vector
    1 + ((sizeof(char) * (IUPAC_FULL_COUNT + 1)) - 1) / sizeof(SparseMotifMapNode) + 1 + ((sizeof(blscounttype) * blsvectorsize) - 1) /  sizeof(SparseMotifMapNode))), // index of first child
range(range) {
    root = new SparseMotifMapNode(startIndexes.first);
    // std::cerr << "motifmap created" << std::endl;
}
void SparseMotifMap::addMotifToMap(const std::string &motif, const int &val) {
    root->addMotifToMap(motif, 0, val, startIndexes, range, blsvectorsize);
}
void SparseMotifMap::recPrintAndDelete(long &unique_count, std::ostream &out) {
    root->recPrintAndDelete("", unique_count, out, startIndexes, range, blsvectorsize);
    // delete root;
    malloc_trim(0); // this gives memory back to OS!
};

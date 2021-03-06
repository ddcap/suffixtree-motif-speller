#ifndef MOTIF_H
#define MOTIF_H


#include <vector>
#include <stack>
#include <iostream>
#include <random>
#include <unordered_set>
#include <bits/stdc++.h>


#define N_BITS 16 // must be at least the maximum number of organisms
typedef unsigned int occurence_bits; // define type here to easily expand number of bits in code!
typedef unsigned short blscounttype;
// typedef unsigned char  blscounttype;


enum IUPAC {
    BASE_A =  0x1, // 0001
    BASE_C =  0x2, // 0010
    IUPAC_M = 0x3, // 0011 A or C
    BASE_G =  0x4, // 0100
    IUPAC_R = 0x5, // 0101 A or G
    IUPAC_S = 0x6, // 0110 G or C
    IUPAC_V = 0x7, // 0111 A or C or G
    BASE_T =  0x8, // 1000
    IUPAC_W = 0x9, // 1001 A or T
    IUPAC_Y = 0xA, // 1010 C or T
    IUPAC_H = 0xB, // 1011 A or C or T
    IUPAC_K = 0xC, // 1100 G or T
    IUPAC_D = 0xD, // 1101 A or G or T
    IUPAC_B = 0xE, // 1110 C or G or T
    IUPAC_N = 0xF  // 1111
};

/**
BLS scsore is calculated like this:
even if leaf the score is added, but 2 leafs in same branch sums the length of leaf but not the branch further
take this part of a tree:
       |------0.3---- A
----0.2|
       |------0.3---- B

i.e. -> 01 = only in A, 10 is only in B and 11 means motif is found in both A and B(we disregard the other leafs for this example)

if only one branch is 0, cause it doesnt connect to others
10 --> score is 0
01 --> score is 0
11 --> score is 0.6!
*/

class Motif {
private:
    static const std::vector<char> complement;

public:
    static std::string getGroupID(const std::string& read);
    static std::string ReverseComplement(const std::string& read);
    static bool isRepresentative(const std::string& read);
    static bool isGroupRepresentative(const std::string& read);
    static std::string getRepresentative(const std::string& read);
    static void writeGroupIDAndMotifInBinary(const std::string& motif, const short &maxlen, std::ostream& out);
    static void writeGroupIDAndMotifInBinary(const long &motifdata, const short &maxlen, std::ostream& out);
    static long getLongRepresentation(const std::string& motif);
    static std::string getStringRepresentation(const long& motif);
    static void writeMotifInBinary(const std::string& motif, const short &maxlen, std::ostream& out);
    static void writeGroupIDAndMotif(const std::string& motif, std::ostream& out);
    static void writeMotif(const std::string& motif, std::ostream& out);
};

// class MotifCollection {
// private:
//     std::unordered_set<std::string> processedMotifs; // average constant time adding + searching!
//
// public:
//     bool checkAndAddElement(std::string motif) {
//         if(processedMotifs.find(motif) == processedMotifs.end()) {
//             processedMotifs.insert(motif);
//             return true;
//         } else {
//             return false;
//         }
//     }
//     size_t size() { return processedMotifs.size(); }
// };



class BLSLinkedListNode {
private:
    float length;
    occurence_bits mask;
    BLSLinkedListNode *next;
    BLSLinkedListNode *child;
    int level;
    std::ostream& write(std::ostream& o) const{
        for(int i = 0; i < level; i++ ) { o << "  "; }
        o << "[" << +mask << "/" << length << "]"; // print this actual node
        if(child != NULL) {
            o << std::endl;
            o << *child;
        }
        if(next != NULL) {
            o << std::endl;
            o << *next;
        }
        return o;
    }
public:
    BLSLinkedListNode(): length(0), mask(0), next(NULL), child(NULL), level(0) {
      for (int i = 0; i < N_BITS; i++) {
        mask |= 1 << i; // so the last bits arent set to 1 for later in popcount etc!;
      }
      // std::cerr << "root node mask: " << +mask << std::endl; // +mask print the actual number not the char!
    }
    BLSLinkedListNode(int level_): length(0), mask(0), next(NULL), child(NULL), level(level_) {}
    BLSLinkedListNode(float length_, occurence_bits mask_, int level_): length(length_), mask(mask_), next(NULL), child(NULL), level(level_) {}
    void setLevel(int level_) { level = level_; }
    BLSLinkedListNode *addNext(int level_) { next = new BLSLinkedListNode(level_); return next;}
    BLSLinkedListNode *addChild(int level_) { child = new BLSLinkedListNode(level_); return child; }
    BLSLinkedListNode *setChild(BLSLinkedListNode *newChild) { child = newChild; return child; }
    void setMask(occurence_bits mask_) {mask = mask_;}
    void setLength(float length_) {length = length_;}
    void addLength(float length_) {length += length_;}

    BLSLinkedListNode *getChild() const { return child; }
    BLSLinkedListNode *getNext() const { return next; }
    occurence_bits getMask() const {return mask; }
    float getLength() const {return length; }

    friend std::ostream& operator<< (std::ostream& o, const BLSLinkedListNode& b) {
        return b.write(o);
    }

    float getScore(const occurence_bits& occurence);
};

class BLSScore {
private:
    std::vector<float> blsThresholds;
    BLSLinkedListNode* root;
    std::vector<float> preparedBLS;
    std::vector<char> preparedBLSVector;

    float calculateBLSScore(const occurence_bits& occurence) const;
    // std::vector<int> calculateBLSVector(const float& bls) const;
    char calculateBLSVector(const float& bls) const;
    void recReadBranch(int recursion, int& leafcount, std::string& newick, BLSLinkedListNode* currentroot, std::vector<std::string> &order_of_species);
    void prepAllCombinations(int used_bits);

public:
    // example: ((BD1G15520:0.2688, OS03G38520:0.2688):0.0538, (SB01G015780:0.086, (ZM01G45380:1.0E-6,ZM05G08300:1.0E-6):0.086):0.2366);
    BLSScore(std::vector<float> blsThresholds_, std::string newick, int species, std::vector<std::string> &order_of_species) : blsThresholds(blsThresholds_){
        root = new BLSLinkedListNode();
        int leafnr = 0;
        recReadBranch(0, leafnr, newick, root, order_of_species);
        prepAllCombinations(species);
    }
    ~BLSScore() {
        // Depth-first traversal of the tree
        std::stack<BLSLinkedListNode *> stack;
        stack.push(root);

        while (!stack.empty()) {
                BLSLinkedListNode* node = stack.top();
                stack.pop();
                if(node->getChild() != NULL) stack.push(node->getChild());
                if(node->getNext() != NULL) stack.push(node->getNext());

                delete node;
        }
    }
    friend std::ostream& operator<< (std::ostream& o, const BLSScore& bls) {
        o << *bls.root << std::endl;
        return o;
    }
    size_t getBLSVectorSize() const { return blsThresholds.size(); }
    float getBLSScore(const occurence_bits& occurence) const;
    const char* getBLSVector(const occurence_bits& occurence) const;
    // const std::vector<int>* getBLSVector(const occurence_bits& occurence) const;
    void writeBLSVectorInBinary(const occurence_bits& occurence, std::ostream& out) const;
    void writeBLSVector(const occurence_bits& occurence, std::ostream& out) const;
    bool greaterThanMinThreshold(const occurence_bits& occurence) const;
    bool greaterThanThreshold(const occurence_bits& occurence, const int& blsThresholdIdx) const;

    blscounttype *createBlsVectorFromByte(const occurence_bits& occurence) const;
    void addByteToBlsVector(blscounttype *v, const occurence_bits& occurence) const;

};

class IupacMask {
private:
    unsigned char mask;
    static const std::vector<std::string> characterLists;

public:
    static const char FILLER = '-';
    static const char DELIMITER = '$';
    static const std::vector<IupacMask> characterToMask;
    static const std::vector<char> representation;
    IupacMask() : mask(0) {}
    IupacMask(const IUPAC mask_ ) : mask(mask_) {}

    unsigned char getMask() const { return mask; }

    bool isDegenerate() {
        return __builtin_popcountll(mask) > 1; // relies on the gcc builtin popcount
    }

    static char getRandomChar(char c);

    const std::string* getCharacters() const;

    char getRepresentation() const;

    friend std::ostream& operator<< (std::ostream& o, const IupacMask& m) {
        o << m.getRepresentation();
        return o;
    }
};


#endif

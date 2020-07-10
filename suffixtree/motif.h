#ifndef MOTIF_H
#define MOTIF_H


#include <vector>
#include <stack>
#include <iostream>
#include <random>
#include <unordered_set>
#include <bits/stdc++.h>


#define N_BITS 8 // must be at least the maximum number of organisms
typedef unsigned short int occurence_bits; // define type here to easily expand number of bits in code!


enum IUPAC {
    BASE_A =  0x1,
    BASE_C =  0x2,
    IUPAC_M = 0x3, // A or C
    BASE_G =  0x4,
    IUPAC_R = 0x5, // A or G
    IUPAC_S = 0x6, // G or C
    IUPAC_V = 0x7, // A or C or G
    BASE_T =  0x8,
    IUPAC_W = 0x9, // A or T
    IUPAC_Y = 0xA, // C or T
    IUPAC_H = 0xB, // A or C or T
    IUPAC_K = 0xC, // G or T
    IUPAC_D = 0xD, // A or G or T
    IUPAC_B = 0xE, // C or G or T
    IUPAC_N = 0xF
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
    static void writeGroupIDAndMotifInBinary(const std::string& motif, std::ostream& out);
    static void writeGroupIDAndMotif(const std::string& motif, std::ostream& out);
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
    void setMask(occurence_bits mask_) {mask = mask_;}
    void setLength(float length_) {length = length_;}

    BLSLinkedListNode *getChild() const { return child; }
    BLSLinkedListNode *getNext() const { return next; }
    occurence_bits getMask() const {return mask; }

    friend std::ostream& operator<< (std::ostream& o, const BLSLinkedListNode& b) {
        return b.write(o);
    }

    float getScore(const occurence_bits& occurence);
};

class BLSScore {
private:
    static const std::vector<float> blsThresholds; // TODO this should be set in constructor or something
    BLSLinkedListNode* root;
    std::vector<float> preparedBLS;
    std::vector<std::vector<int> > preparedBLSVector;

    float calculateBLSScore(const occurence_bits& occurence) const;
    std::vector<int> calculateBLSVector(const float& bls) const;
    void recReadBranch(int recursion, int& leafcount, std::string& newick, BLSLinkedListNode* currentroot);
    void prepAllCombinations();

public:
    // example: ((BD1G15520:0.2688, OS03G38520:0.2688):0.0538, (SB01G015780:0.086, (ZM01G45380:1.0E-6,ZM05G08300:1.0E-6):0.086):0.2366);
    BLSScore(std::string newick) {
        root = new BLSLinkedListNode();
        int leafnr = 0;
        recReadBranch(0, leafnr, newick, root);
        prepAllCombinations();
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
    float getBLSScore(const occurence_bits& occurence) const;
    const std::vector<int>* getBLSVector(const occurence_bits& occurence) const;
    void writeBLSVectorInBinary(const occurence_bits& occurence, std::ostream& out) const;
    void writeBLSVector(const occurence_bits& occurence, std::ostream& out) const;
    char readBLSVectorInBinary(std::istream& in) const;
    bool greaterThanMinThreshold(const occurence_bits& occurence) const;
};

class IupacMask {
private:
    unsigned char mask;
    static const std::vector<std::string> characterLists;
    static const std::vector<char> representation;

public:
    static const char FILLER = '-';
    static const char DELIMITER = '$';
    static const std::vector<IupacMask> characterToMask;
    IupacMask() : mask(0) {}
    IupacMask(const IUPAC mask_ ) : mask(mask_) {}

    unsigned char getMask() const { return mask; }

    bool isDegenerate() {
        return __builtin_popcountll(mask) > 1; // relies on the gcc builtin popcount
    }

    const std::string* getCharacters() const;

    char getRepresentation() const;

    friend std::ostream& operator<< (std::ostream& o, const IupacMask& m) {
        o << m.getRepresentation();
        return o;
    }
};


#endif

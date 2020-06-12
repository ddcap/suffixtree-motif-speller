#ifndef MOTIF_H
#define MOTIF_H


#include <vector>
#include <stack>
#include <iostream>
#include <bitset>
#include <random>


#define N_BITS 4

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

template<unsigned char N>
class BLSLinkedListNode {
private:
    float length;
    std::bitset<N> mask;
    BLSLinkedListNode<N> *next;
    BLSLinkedListNode<N> *child;
    int level;
    std::ostream& write(std::ostream& o) const;
public:
    BLSLinkedListNode<N>(): length(0), mask(0), next(NULL), child(NULL), level(0) {
        mask.flip(); // root node (all 1!)
    }
    BLSLinkedListNode<N>(int level_): length(0), mask(0), next(NULL), child(NULL), level(level_) {}
    ~BLSLinkedListNode<N>() {
        if (child != NULL) { delete child; }
        if (next != NULL) { delete next; }
    }
    BLSLinkedListNode<N>(float length_, std::bitset<N> mask_, int level_): length(length_), mask(mask_), next(NULL), child(NULL), level(level_) {}
    void setLevel(int level_) { level = level_; }
    BLSLinkedListNode<N> *addNext(int level_) { next = new BLSLinkedListNode<N>(level_); return next;}
    BLSLinkedListNode<N> *addChild(int level_) { child = new BLSLinkedListNode<N>(level_); return child; }
    void setMask(std::bitset<N> mask_) {mask = mask_;}
    void setLength(float length_) {length = length_;}

    BLSLinkedListNode<N> *getChild() { return child; }
    BLSLinkedListNode<N> *getNext() { return next; }
    std::bitset<N> getMask() {return mask; }

    friend std::ostream& operator<< (std::ostream& o, const BLSLinkedListNode<N>& b) {
        return b.write(o);
    }

    float getScore(const std::bitset<N>& occurence);
};

class BLSScore {
private:
    BLSLinkedListNode<N_BITS>* root;

    void recReadBranch(int recursion, int& leafcount, std::string& newick, BLSLinkedListNode<N_BITS>* currentroot);

public:
    // example: ((BD1G15520:0.2688, OS03G38520:0.2688):0.0538, (SB01G015780:0.086, (ZM01G45380:1.0E-6,ZM05G08300:1.0E-6):0.086):0.2366);
    BLSScore(std::string newick) {
        root = new BLSLinkedListNode<N_BITS>();
        int leafnr = 0;
        recReadBranch(0, leafnr, newick, root);
    }
    ~BLSScore() {
        delete root;
    }
    friend std::ostream& operator<< (std::ostream& o, const BLSScore& bls) {
        o << *bls.root << std::endl;
        return o;
    }
    float getBLSScore(std::bitset<N_BITS> occurence) const;
};

class IupacMask {
private:
    unsigned char mask;

public:
    IupacMask(const IUPAC mask_ ) : mask(mask_) {

    }

    bool isDegenerate() {
        return __builtin_popcountll(mask) > 1;
    }

    void getCharacters(std::vector<char>& chars);

    char getRepresentation() const;

    friend std::ostream& operator<< (std::ostream& o, const IupacMask& m) {
        o << m.getRepresentation();
        return o;
    }
};


#endif

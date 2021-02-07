/***************************************************************************
 *   Copyright (C) 2018 Jan Fostier (jan.fostier@ugent.be)                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef SUFFIXTREE_H
#define SUFFIXTREE_H

#include <tuple>
#include <limits>
#include <string>
#include <vector>
#include <unordered_set>
#include <cassert>
#include "motif.h"
#include "motifmap.h"

// ============================================================================
// (TYPE) DEFINITIONS AND PROTOTYPES
// ============================================================================

typedef uint32_t length_t;
//typedef unsigned short blscounttype;
typedef unsigned char blscounttype;
typedef SparseMotifMap MyMotifMap;

// <position in T, position in Q, length of MEM>
typedef std::tuple<size_t, size_t, size_t> MEMOcc;

enum Alphabet { EXACT = 0x0, EXACTANDN = 0x1, TWOFOLDSANDN = 0x2, ALL = 0x3 };

struct MotifPosition {
  int family;
  int gene;
  int reverseComplement;
  int position;
};

#define MAX_ASCII_CHAR 85
#define MAX_CHAR 7
// #define MAX_CHAR 256
// 7 characters in string -> acgt n $ and ' '(space) others will give index -1 in array which will crash the program.

// ============================================================================
// CLASS SUFFIX TREE NODE
// ============================================================================

// A suffix tree node contains links to its parent and children, a suffix link,
// the depth of the node and a suffix index in case it is a leaf node.
// For non-leaf nodes, the suffix index contains the length_t::max() value.

// A suffix tree node also contains a range [beginIdx, endIdx[ in T of its
// parent edge. The range encodes the characters implied on the edge.
class STNode {

private:
        // edge properties (between parent and current node)
        length_t beginIdx;              // begin index in T of parent edge
        length_t endIdx;                // end index in T of parent edge

        // node properties
        STNode* parent;                 // pointer to the parent (NULL for root)
        STNode* child[MAX_CHAR];        // pointers to children
        STNode* suffixLink;             // points to suffix link node
        length_t depth;                 // depth of current node
        length_t suffixIdx;             // suffix index (only for leaf nodes)
        occurence_bits occurence;
        static const std::vector<char> Alphabet;
        // static const std::vector<short> charToIndex;
        static const short charToIndex[MAX_ASCII_CHAR];

public:
        /**
         * Constructor
         * @param begin Begin index in T
         * @param end End index in T
         */
        STNode(length_t begin, length_t end) : beginIdx(begin), endIdx(end), occurence(0) {
                parent = NULL;
                for (int i = 0; i < MAX_CHAR; i++)
                        child[i] = NULL;
                suffixLink = NULL;
                depth = end - begin;
                suffixIdx = std::numeric_limits<length_t>::max();
        }


        /**
         * Sets the bit for string number 'occurenceBit' to true in a GST
         * @return occurenceBit is the number of the current string this suffix belongs to
         */
        void setOccurenceBitForGST(unsigned char occurenceBit) { // used in leaf
            occurence |= 1 << (occurenceBit ); // divide by 2 to take into account Reverse complement!
            if(parent != NULL && !(parent->occurence &  (1 << (occurenceBit )) ) ) { // if parent and parent hasnt got this occurence set it!
              parent->setOccurenceBitForGST(occurenceBit);
            }
        }
        void setOccurence(occurence_bits occurence_) {
            occurence = occurence_;
        }
        occurence_bits getOccurence() {
          return occurence;
        }

        /**
         * Get the begin index in T of parent edge
         * @return Begin index
         */
        length_t begin() const {
                return beginIdx;
        }

        /**
         * Set the begin index in T of parent edge
         * @param target target value
         */
        void setBegin(length_t target) {
                beginIdx = target;
        }

        /**
         * Get the end index in T of parent edge
         * @return End index
         */
        length_t end() const {
                return endIdx;
        }

        /**
         * Set the begin index in T of parent edge
         * @param target target value
         */
        void setEnd(length_t target) {
                endIdx = target;
        }

        /**
         * Get the length of the edge between parent en this node
         * @return Edge length
         */
        length_t getEdgeLength() const {
                return endIdx - beginIdx;
        }

        /**
         * Get the parent node
         * @return The parent node
         */
        STNode* getParent() const {
                return parent;
        }

        /**
         * Get a pointer to child node for which the edge starts with c
         * @param c Character c
         */
        STNode* getChild(char c) const {
                return child[charToIndex[static_cast<unsigned char>(c)]];
                // return child[static_cast<unsigned char>(c)];
        }
        STNode* getChildNumber(unsigned char i) const {
                return child[i];
        }

        /**
         * Set a (pre-allocated) child for this node
         * @param c First character on the edge to the child
         * @param chdToAdd Pointer to the child
         */
        void setChild(char c, STNode* chdToAdd) {
                chdToAdd->parent = this;
                chdToAdd->depth = depth + chdToAdd->getEdgeLength();
                child[charToIndex[static_cast<unsigned char>(c)]] = chdToAdd;
                // child[static_cast<unsigned char>(c)] = chdToAdd;
        }

        /**
         * Set the suffix link (for internal nodes only)
         * @param target Target value
         */
        void setSuffixLink(STNode* target) {
                suffixLink = target;
        }

        /**
         * Get the suffix link from this node (only for internal nodes)
         * @return suffix link
         */
        STNode* getSuffixLink() const {
                return suffixLink;
        }

        /**
         * Set the suffix index (for leaves only)
         * @param target Target value
         */
        void setSuffixIdx(length_t target) {
                suffixIdx = target;
        }

        /**
         * Get the suffix index in T (for leaves only)
         * @return suffix index
         */
        length_t getSuffixIdx() const {
                return suffixIdx;
        }

        /**
         * Set the depth of the node
         * @param target Target value
         */
        void setDepth(length_t target) {
                depth = target;
        }

        /**
         * Get the depth of the node
         * @return Depth of node
         */
        length_t getDepth() const {
                return depth;
        }

        /**
         * Check whether the node is a leaf
         * @return true or false
         */
        bool isLeaf() const {
                return suffixIdx != std::numeric_limits<length_t>::max();
        }
};

// ============================================================================
// CLASS SUFFIX TREE POSITION
// ============================================================================

// A position in a suffix tree is a pointer to some node of the (implied)
// uncompressed suffix trie. In other words, a position could point to a leaf
// or an internal node in the suffix tree, or to a position halfway an edge.
// In that case, the position is specified by a pointer to the child node of
// that edge and an offset along the edge.

class STPosition {

private:
        STNode* node;                   // pointer to child node
        length_t offset;                // offset > 0 along the edge

public:
        /**
         * Constructor for position that points to a node itself
         * @param node Pointer to node in the ST
         */
        STPosition(STNode* node) : node(node), offset(node->getEdgeLength()) {}

        /**
         * Constructor for position that points to an edge
         * @param node Pointer to the child node
         * @param offset Offset along edge from parent to child
         */
        STPosition(STNode* node, length_t offset) : node(node), offset(offset) {
                assert(offset <= node->getEdgeLength());
        }
        void set(STNode* node_)  {
                node = node_;
                offset = node->getEdgeLength();
                assert(offset <= node->getEdgeLength());
        }
        void set(STNode* node_, length_t offset_)  {
                node = node_;
                offset = offset_;
                assert(offset <= node->getEdgeLength());
        }

        occurence_bits getOccurence() {
          return node->getOccurence();
        }
        /**
         * Check whether position points to a node (true) or edge (false)
         * @return True or false
         */
        bool atNode() const {
                return offset == node->getEdgeLength();
        }

        size_t getPositionInText() const {
          return node->getParent() == NULL ? -1 : node->begin() + offset; // -1 for root node, else position in string
        }
        /**
         * Get the depth of the position
         * @return Depth of the position
         */
        length_t getDepth() const {
                // do not use node->getParent()->depth as root has no parent
                return node->getDepth() - node->getEdgeLength() + offset;
        }

        /**
         * Operator== overloading
         * @param rhs Right hand side
         * @return true or false
         */
        bool operator==(const STPosition& rhs) const {
                if (node != rhs.node)
                        return false;
                if (offset != rhs.offset)
                        return false;
                return true;
        }

        /**
         * Operator != overloading
         * @param rhs Right hand side
         * @return true or false
         */
        bool operator!=(const STPosition& rhs) {
                return !(*this == rhs);
        }

        friend class SuffixTree;
};

// per character you add a list of positions
struct STPositionVector {
public:
  std::vector<STPosition> list;
  size_t validPositions;
  STPositionVector(int maxDenegeracy): validPositions(0) {
    list.reserve(pow(4, maxDenegeracy)); // max number of positions is 4^maxDegen, ACGT = 4 on maxDegen positions!
  }
  void reset() {validPositions = 0;}
  void addSTPosition(STNode* node, length_t offset) {
    list[validPositions].set(node, offset);
    validPositions++;
  }
  void addSTPosition(STNode* node) {
    list[validPositions].set(node);
    validPositions++;
  }
  bool empty() {
    return validPositions == 0;
  }
};

// the list for every character
struct STPositionsPerLetter {
public:
  std::vector<STPositionVector> list; // PUBLIC SO IT ISNT COPIED ALL THE TIME!
  STPositionsPerLetter(int motifSizeUpperBound, int maxDenegeracy) {
    for(int i= 0; i < motifSizeUpperBound; i++) { // one more for last iteration!
      list.push_back(STPositionVector(maxDenegeracy));
    }
  }
  void reset() {
    for(size_t i= 0; i < list.size(); i++) {
      list[i].reset();
    }
  }
};

// ============================================================================
// CLASS SUFFIX TREE
// ============================================================================

class SuffixTree;
typedef void (SuffixTree::*printMotifPtr)(const short& maxlen, const std::string& currentMotif, const BLSScore& bls, const occurence_bits& occurence, std::ostream& out);

class SuffixTree {

private:
        // --------------------------------------------------------------------
        // ROUTINES TO MANIPULATE SUFFIX TREE POSITIONS
        // --------------------------------------------------------------------

        /**
         * Given a position, try and advance by matching a single character
         * @param pos Position to advance (input / output)
         * @param c Character c to match
         * @return true if position was advanced, false otherwise
         */
        bool advancePos(STPosition& pos, char c) const;

        /**
         * Given a position, try and advance by matching a pattern
         * This procedure advances a much as possible, even if P cannot be
         * entirely matched
         * @param pos Position to advance (input / output)
         * @param P Pattern to match str[begin, end[
         * @param begin First position in str to match
         * @param end End position in str to match
         * @return true if pattern is fully matched, false otherwise
         */
        bool advancePos(STPosition& pos, const std::string& P,
                        size_t begin, size_t end) const;

        /**
         * Faster matching using the skip/count trick (patterns *must* exist!)
         * @param pos Position to advance (input / output)
         * @param P Pattern to match str[begin, end[
         * @param begin First position in str to match
         * @param end End position in str to match
         */
        void advancePosSkipCount(STPosition& pos, const std::string& P,
                                 size_t begin, size_t end) const;

        /**
         * Given a suffix tree position, get the position by following the SL
         * @param pos Input position
         * @return Position after following the suffix link
         */
        STPosition followSuffixLink(const STPosition& pos) const;

        /**
         * Get the implied string from a suffix tree position
         * @param pos Suffix tree position
         * @return Implied string from root to pos
         */
        std::string posToStr(const STPosition& pos) const;

        // --------------------------------------------------------------------
        // ROUTINES FOR PATTERN MATCHING
        // --------------------------------------------------------------------

        /**
         * Get all occurrences under a given position
         * @param node Pointer to a node
         * @param occ Vector to store the occurrences (output)
         */
        void getOccurrences(const STPosition& pos, std::vector<size_t>& occ) const;

        /**
         * Find the Maximal Exact Matches between T and P
         * @param Q Query pattern
         * @param i Current suffix of Q under consideration
         * @param minSize Minimum size of the MEM
         * @param pos Position in ST that points to RMEM of suf_i(Q)
         * @param occ MEM occurrences (output)
         */
        void reportMEM(const std::string& Q, size_t i, size_t minSize,
                       const STPosition& pos, std::vector<MEMOcc>& occ);

        // --------------------------------------------------------------------
        // ROUTINES TO CONSTRUCT/MANIPULATE SUFFIX TREE
        // --------------------------------------------------------------------

        /**
         * Given a suffix tree position, split the edge
         * @param pos Suffix tree position
         * @return Position of the edge has been split
         */
        STPosition splitEdge(const STPosition& pos);

        /**
         * Add a leaf to the suffix tree
         * @param pos Suffix tree position
         * @param suffixIndex
         */
        void addLeaf(const STPosition& pos, length_t suffixIndex);
        void addLeaf(const STPosition& pos, length_t suffixIndex, unsigned char currentbit);

        /**
         * Step 1 of the Maass algorithm for the computation of suffix links
         * @param node Current node
         * @param A Vector that maps causing leaf idx to the source of a SL
         * @return The minimum suffix index for the subtree under node
         */
        length_t recComputeSLPhase1(STNode* node, std::vector<STNode*>& A);

        /**
         * Step 2 of the Maass algorithm for the computation of suffix links
         * @param node Current node
         * @param A Vector that maps causing leaf idx to the source of a SL
         * @param B Vector that maps depth to node on the current path
         */
        void recComputeSLPhase2(STNode* node, const std::vector<STNode*>& A,
                                std::vector<STNode*>& B);

        /**
         * Compute suffix links in linear time using the Maass algorithm
         */
        void computeSuffixLinks();

        /**
         * Construct the suffix tree using a naive O(n^2) time algorithm
         */
        void constructNaive();

        /**
         * Construct the suffix tree using Ukonen's algorithm in O(n) time
         */
        void constructUkonen();

        // --------------------------------------------------------------------
        // ROUTINES FOR I/O
        // --------------------------------------------------------------------

        /**
         * Write suffix tree to the output stream in a formatted way
         * @param o Output stream
         * @return Output stream
         */
        std::ostream& write(std::ostream& o) const;

        // --------------------------------------------------------------------
        const std::string T;            // text to index
        const std::string name;            // text to index
        STNode* root;                   // pointer to the root node
        static const std::vector<IupacMask> exactAlphabet;
        static const std::vector<IupacMask> exactAndNAlphabet;
        static const std::vector<IupacMask> exactTwofoldAndNAlphabet;
        static const std::vector<IupacMask> exactAndAllDegenerateAlphabet;
        const std::vector<IupacMask> *alphabet;
        int reverseComplementFactor = 1;
        int motifCount;
        size_t iteratorCount;
        size_t node_count = 0;
        std::vector<size_t> stringStartPositions; // indicates where new strings start
        std::vector<std::string> gene_names; // identify gene names
        std::vector<size_t> next_gene_locations; // identify genes
        std::vector<size_t> order_of_species_mapping; // map species to correct index in the bls tree
        MyMotifMap *motifmap = NULL;
        // --------------------------------------------------------------------

        void recPrintMotifs(const std::pair<short, short>& l,
          const int& maxDegenerateLetters, const BLSScore& bls,
          STPositionsPerLetter& matchingNodes, const std::string& prefix,
          int curDegenerateLetters, std::ostream& out);
        // this next one also returns all positions the current Motif matches
        void recPrintMotifsWithPositions(const std::pair<short, short>& l,
          const int& maxDegenerateLetters, const BLSScore& bls,
          STPositionsPerLetter& matchingNodes, std::vector<std::pair<int, int>>& stringPositions, const std::string& prefix,
          int curDegenerateLetters, std::ostream& out);

        void getLeafPositions(std::vector<std::pair<int, int>>& positions, const std::vector<STPosition>& nodePositions, const size_t size) const;
        void getPositionsStartingWithDelimiter(std::vector<std::pair<int, int>>& positions, const std::vector<STPosition>& nodePositions, const size_t size) const;

        void advanceIupacCharacter(const IupacMask& mask, const int& characterPos, STPositionsPerLetter& matchingNodes, occurence_bits& occurence) const;
        void advanceExactCharacter(const IupacMask& mask, const int& characterPos, STPositionsPerLetter& matchingNodes, occurence_bits& occurence) const;
        void getBestOccurence(std::vector<std::pair<int, int>>& positions, const BLSScore& bls, occurence_bits& occurence);

        void printMotifBinary(const short& maxlen, const std::string& currentMotif, const BLSScore& bls, const occurence_bits& occurence, std::ostream& out);
        void printMotifString(const short& maxlen, const std::string& currentMotif, const BLSScore& bls, const occurence_bits& occurence, std::ostream& out);

        void addMotifToMap(const short& maxlen, const std::string& currentMotif, const BLSScore& bls, const occurence_bits& occurence);

        printMotifPtr printMotif = &SuffixTree::printMotifBinary;
        // printMotifPtr printMotif = &SuffixTree::printMotifString;

        static bool positionPairSort(const std::pair<int,int> &a,  const std::pair<int,int> &b)
        {
            return a.second < b.second;
        }

public:
        /**
         * Constructor
         * @param T Text to be indexed
         */
        // SuffixTree(const std::string& T) : SuffixTree(T, false) {}
        // SuffixTree(const std::string& T, bool hasReverseComplement);
        SuffixTree(const std::string& T, const std::string& name, bool hasReverseComplement, std::vector<size_t> stringStartPositions_, std::vector<std::string> gene_names_,
          std::vector<size_t> next_gene_locations_, std::vector<size_t> order_of_species_mapping_, MyMotifMap *motifmap_);
        //

        /**
         * Destructor
         */
        ~SuffixTree();

        /**
         * Find occurrences of P in the suffix tree in O(m + #occ) time
         * @param P Pattern P to match (length m)
         * @param occ Start positions of the occurrences in T (output)
         */
        void matchPattern(const std::string& P, std::vector<size_t>& occ);
        std::vector<std::pair<int, int>> matchIupacPattern(const std::string& P, const BLSScore& bls, int maxDegenerateLetters, occurence_bits& occurence);
        std::vector<std::pair<int, int>> matchIupacPatternWithPositions(const std::string& P, const BLSScore& bls, int maxDegenerateLetters, occurence_bits& occurence);

        int matchIupacPatterns(std::istream& in, std::ostream& out, const BLSScore& bls, const int &maxDegenerateLetters, const short& maxlen);
        // TODO add same but with positions
        void matchPattern(const std::string& P, BLSScore& bls);
        void printMotifPositions(std::ostream& out, const std::string &motif, std::vector<std::pair<int, int>> positions, const float blsScore);
        void getLeafPositionsAndPrint(const std::vector<STPosition>& matchingNodes, const size_t size,
          std::ostream& out, const std::string &motif, const float blsScore) const;
        size_t getMotifsIteratedCount() { return iteratorCount; }


        /**
         * Find the Maximal Exact Matches between T and P in O(m + #RMEMs) time
         * @param Q Query pattern Q (length m)
         * @param minSize Minimum size of the MEM
         * @param occ MEM occurrences (output)
         */
        void findMEM(const std::string& Q, size_t minSize,
                     std::vector<MEMOcc>& occ);

        /**
         * Operator<< overloading
         * @param o Output stream
         * @param ST Suffix tree
         * @return Output stream with ST information
         */
        friend std::ostream& operator<< (std::ostream& o, const SuffixTree& t) {
                return t.write(o);
        }

        /**
        * Find all motifs in the Suffix tree of length l
        * @Param l Length of motifs to find
        * @Param motifs STPositions of different motifs in T (output)
        */
        int printMotifs(const std::pair<short, short>& l, const Alphabet alphabet,
          const int& maxDegenerateLetters, const BLSScore& bls,
          std::ostream& out, bool isAlignmentBased);

};

#endif

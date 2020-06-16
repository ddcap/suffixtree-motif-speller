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

#include <iostream>
#include <cassert>
#include <stack>
#include "suffixtree.h"
#include "motif.h"

using namespace std;

// ============================================================================
// SUFFIX TREE (PRIVATE FUNCTIONS)
// ============================================================================

// ----------------------------------------------------------------------------
// ROUTINES TO MANIPULATE SUFFIX TREE POSITIONS
// ----------------------------------------------------------------------------

bool SuffixTree::advancePos(STPosition& pos, char c) const
{
        // case a) we are at a node: find a new child
        if (pos.atNode()) {
                STNode* chd = pos.node->getChild(c);
                if (chd == NULL)        // no such edge, get out
                        return false;
                pos.node = chd;
                pos.offset = 1;
                return true;
        }

        // case b) we are at an edge: try and match next character along edge
        if (T[pos.node->begin() + pos.offset] != c)
                return false;

        pos.offset++;
        return true;
}

bool SuffixTree::advancePos(STPosition& pos, const string& P,
                            size_t begin, size_t end) const
{
        for (auto itP = P.begin() + begin; itP < P.begin() + end; itP++)
                if (!advancePos(pos, *itP))
                        return false;

        return true;
}

void SuffixTree::advancePosSkipCount(STPosition& pos, const std::string& P,
                                     size_t begin, size_t end) const
{
        assert(begin <= end);

        // progress as much as possible on the current edge
        length_t adv = min<length_t>(end-begin, pos.node->getEdgeLength()-pos.offset);
        pos.offset += adv;
        begin += adv;

        while (begin < end) {
                // move to the correct child
                pos.node = pos.node->getChild(P[begin]);

                // progress as much as possible on the current edge
                adv = min<length_t>(end-begin, pos.node->getEdgeLength());
                pos.offset = adv;
                begin += adv;
        }
}

STPosition SuffixTree::followSuffixLink(const STPosition& pos) const
{
        // for root, do nothing
        if (pos.node == root)
                return pos;

        STNode* parent = pos.node->getParent();
        if (parent == root) {   // special case for edges originating from the root
                STPosition newPos(root);
                advancePosSkipCount(newPos, T, pos.node->begin()+1, pos.node->begin()+pos.offset);
                return newPos;
        } else {                // generic case (parent is an internal node)
                STPosition newPos(parent->getSuffixLink());
                advancePosSkipCount(newPos, T, pos.node->begin(), pos.node->begin()+pos.offset);
                return newPos;
        }
}

string SuffixTree::posToStr(const STPosition& pos) const
{
        string str;
        str.reserve(pos.getDepth());

        // create a stack with nodes from position (bottom) to root (top)
        stack<STNode*> path;
        path.push(pos.node);
        while (path.top()->getParent() != NULL)
                path.push(path.top()->getParent());

        path.pop();     // remove root (empty string);

        // traverse the stack from root (top) to position (bottom)
        while(!path.empty()) {
                STNode* node = path.top();
                path.pop();

                size_t len = path.empty() ? pos.offset : node->getEdgeLength();
                str.append(T.substr(node->begin(), len));
        }

        return str;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR PATTERN MATCHING
// ----------------------------------------------------------------------------

void SuffixTree::getOccurrences(const STPosition& pos, vector<size_t>& occ) const
{
        // if so, find all the leaves under "pos"
        stack<STNode*> stack;
        stack.push(pos.node);

        while (!stack.empty()) {
                STNode* node = stack.top();
                stack.pop();

                if (node->isLeaf())
                        occ.push_back(node->getSuffixIdx());
                else
                        for (int i = 0; i < MAX_CHAR; i++)
                                if (node->getChild(i) != NULL)
                                        stack.push(node->getChild(i));
        }
}

void SuffixTree::reportMEM(const string& Q, size_t j, size_t minSize,
                           const STPosition& pos, vector<MEMOcc>& occ)
{
        vector<size_t> this_occ;
        getOccurrences(pos.node, this_occ);
        for (auto i : this_occ)
                if (j == 0 || i == 0 || Q[j-1] != T[i-1]) // left-maximal?
                        occ.push_back(MEMOcc(i, j, pos.getDepth()));

        STNode* node = pos.node->getParent();
        STNode* last = pos.node;

        while (node->getDepth() >= minSize) {
                for (size_t i = 0; i < MAX_CHAR; i++) {
                        STNode* chd = node->getChild(i);
                        if (chd == NULL || chd == last)
                                continue;

                        vector<size_t> this_occ;
                        getOccurrences(chd, this_occ);

                        for (auto i : this_occ)
                                if (i == 0 || j == 0 || Q[j-1] != T[i-1]) // left-maximal?
                                        occ.push_back(MEMOcc(i, j, node->getDepth()));
                }

                last = node;
                node = node->getParent();
        }
}

// ----------------------------------------------------------------------------
// ROUTINES TO MANIPULATE SUFFIX TREE
// ----------------------------------------------------------------------------

STPosition SuffixTree::splitEdge(const STPosition& pos)
{
        // check whether we are really in the middle of an edge
        if (pos.atNode())
                return pos;

        STNode *chd = pos.node;
        STNode *par = chd->getParent();
        STNode *mid = new STNode(chd->begin(), chd->begin() + pos.offset);
        mid->setOccurence(chd->getOccurence());
        chd->setBegin(mid->end());


        // set the correct pointers
        par->setChild(T[mid->begin()], mid);
        mid->setChild(T[chd->begin()], chd);

        return STPosition(mid);
}

void SuffixTree::addLeaf(const STPosition& pos, length_t suffixIndex)
{
        STNode *leaf = new STNode(suffixIndex + pos.getDepth(), T.size());
        leaf->setSuffixIdx(suffixIndex);
        pos.node->setChild(T[suffixIndex + pos.getDepth()], leaf);
        // std::cerr << "added new leaf:" << T.substr(suffixIndex + pos.getDepth()) << " ";
}
void SuffixTree::addLeaf(const STPosition& pos, length_t suffixIndex, unsigned char currentbit)
{
        addLeaf(pos, suffixIndex);
        pos.node->getChild(T[suffixIndex + pos.getDepth()])->setOccurenceBitForGST(currentbit, reverseComplementFactor); // factor 2 means reversecomplement is added in the reference string!
}

length_t SuffixTree::recComputeSLPhase1(STNode* node, vector<STNode*>& A)
{
        // for leaves, simply return the suffix index
        if (node->isLeaf())
                return node->getSuffixIdx();

        length_t min1 = T.length() + 1; // the smallest suffix idx at this node
        length_t min2 = T.length() + 1; // the second smallest suffix idx

        for (size_t i = 0; i < MAX_CHAR; i++) {
                if (node->getChild(i) == NULL)
                        continue;

                length_t m = recComputeSLPhase1(node->getChild(i), A);
                if (m < min1) {
                        min2 = min1;
                        min1 = m;
                } else if (m < min2) {
                        min2 = m;
                }
        }

        // indicates that there is a suffix link from "node" to the internal
        // node that is caused by leaf with suffix index min2 + 1
        A[min2 + 1] = node;

        return min1;
}

void SuffixTree::recComputeSLPhase2(STNode* node, const vector<STNode*>& A,
                                    vector<STNode*>& B)
{
        if (node->isLeaf()) {
                length_t sufIdx = node->getSuffixIdx();
                if (A[sufIdx] != NULL && A[sufIdx] != root) {
                        length_t d = A[sufIdx]->getDepth();
                        A[sufIdx]->setSuffixLink(B[d-1]);
                }
        } else {
                length_t d = node->getDepth();
                B[d] = node;

                for (size_t i = 0; i < MAX_CHAR; i++)
                        if (node->getChild(i) != NULL)
                                recComputeSLPhase2(node->getChild(i), A, B);
        }
}

void SuffixTree::computeSuffixLinks()
{
        // see paper by Moritz G. Maass for definintions of cause() and branch()

        // compute A such that A[cause(node) + 1] = node
        vector<STNode*> A(T.size(), NULL);
        recComputeSLPhase1(root, A);

        // compute B such that B[depth] = branch(node, depth)
        vector<STNode*> B(T.size(), NULL);
        recComputeSLPhase2(root, A, B);
}

// ----------------------------------------------------------------------------
// ROUTINES FOR I/O
// ----------------------------------------------------------------------------

std::ostream& SuffixTree::write(std::ostream& o) const
{
        stack<pair<int, STNode*> > stack;
        stack.push(make_pair(0, root));

        while (!stack.empty()) {
                int depth = stack.top().first;
                STNode* node = stack.top().second;
                stack.pop();

                for (int i = 0; i < depth; i++)
                        cout << "  ";
                string nodeDescr = (node->isLeaf())? "LEAF " : "INTL ";
                if (node == root)
                        nodeDescr = "ROOT ";

                o << nodeDescr << "\""
                  << T.substr(node->begin(), node->getEdgeLength()) << "\""
                  << ", depth=" << node->getDepth() << ", occ: " << node->getOccurence();

                if (!node->isLeaf() && node != root && node->getSuffixLink() != NULL)
                        o << " (SL: \"" << posToStr(node->getSuffixLink()) << "\")";
                if (node->isLeaf())
                        o << " (" << node->getSuffixIdx() << ")";
                o << "\n";

                for (int i = MAX_CHAR-1; i >= 0; i--) {
                        char c = (char)i;
                        if (node->getChild(c) != NULL)
                                stack.push(make_pair(depth+1, node->getChild(c)));
                }
        }

        return o;
}

void SuffixTree::constructNaive()
{
        // create a node with an empty range (it has no parent)
        root = new STNode(0, 0);

        for (size_t i = 0; i < T.size(); i++) {
                // find position in ST that maximally matches the prefix of suf_i(T)
                STPosition pos(root);

                // if suf_i(T) is completely matched, T is not prefix-free
                if (advancePos(pos, T, i, T.size()))
                        break;
                        //throw runtime_error("Text T is not prefix-free");

                // middle of an edge? -> first split edge first
                if (!pos.atNode())
                        pos = splitEdge(pos);

                // add a new leaf to the suffix tree
                addLeaf(pos, i);
        }

        // Maass algorithm to compute suffix links in O(n) time
        computeSuffixLinks();
}

/**
I add occurence tables in every leaf node, while keeping the parents satisfied (if bit is not set in parent set in parent, until parent == null or bit is set)
When we split a branch, I copy the bitset of the child (can be a leaf or a node) to the mid node (parent should already be fine since this isnt a new node)
And then add the new bit to the addleaf as we do normally, recursing through parents.
*/
void SuffixTree::constructUkonen()
{
        // create root node with an empty range (it has no parent)
        root = new STNode(0, 0);

        // algorithm invariant: pos points to T[i:j-1[
        STPosition pos(root);
        unsigned char currentleafbit = 0; // supports up to 255 strings!

        // in phase j, build implicit suffix tree for prefix T[0:j[
        for (size_t j = 1, numLeaves = 0; j <= T.size(); j++) {
                // std::cout << "current j " << j << " is at " << T[j] << std::endl;
                STNode *prevInternal = NULL;
                // i starts from numleaves -> keep the current bit to set in a variable!

                // skip 'numleafs' times rule 1 (extension of leaf)
                unsigned char currentbit = currentleafbit;
                for (size_t i = numLeaves; i < j; i++) {
                        // std::cout << "current i " << i << " is at " << T[i] << std::endl;
                        // note that pos will always point to T[i:j-1[ at this point
                        // TODO add gst occurency table
                        // if t[i] == # then set bit of currentbit + increment currentbit
                        // add bit to existing nodes


                        // add a SL from the previously created internal node
                        if (prevInternal != NULL && pos.atNode()) {
                                prevInternal->setSuffixLink(pos.node);
                                prevInternal = NULL;
                        }

                        // try and add character T[j-1]
                        if (advancePos(pos, T[j-1])) {
                            if(numLeaves != i && T[i] == '#') { // only increment one if i == numleaves
                                // std::cout << "incrementing currentbit at " << T[i] << std::endl;
                                currentbit++;
                            } // increment currentbit for next iteration in for loop.
                            // std::cerr << "break at " << i << std::endl;
                            break; // rule 3 : do nothing + show stopper

                        }

                        // rule 2 : create internal node (optionally) and leaf
                        if (!pos.atNode()) {
                                pos = splitEdge(pos); // pos points to new node // TODO add currentbit variable
                                if (prevInternal != NULL)
                                        prevInternal->setSuffixLink(pos.node);
                                prevInternal = pos.node;
                        }

                        addLeaf(pos, i, currentbit); // TODO add currentbit variable
                        if(T[numLeaves] == '#') {
                            // std::cout << "incrementing currentleafbit at " << T[numLeaves] << std::endl;
                            currentleafbit++;
                        } // if previous character where the suffix starts is # then a new string has started so increment currentleafbit
                        if(T[i] == '#') {
                            // std::cout << "incrementing currentbit at " << T[i] << std::endl;
                            currentbit++;
                        } // increment currentbit for next iteration in for loop.

                        // extension is complete: follow suffix link
                        pos = followSuffixLink(pos);

                        numLeaves++;
                        // std::cerr << *this << std::endl;
                }
        }
}

// Routines to explore SuffixTree

/**

Recurse into the tree with degenerate letters, meaning we need a vector of positions we are currently at!
with 3 N's the most number of positions is 4*4*4= 64 positions to check.
To check validity of the motif with regard to the BLS score, we do the or operation on the mask of every position and then check this occurence bitset in the bls score function.
If score is equal to or more than the threshold then we keep the motif if its less we should stop the recursion!
*/
void SuffixTree::recPrintMotifs(const std::pair<short, short>& l,
    const int& maxDegenerateLetters, const BLSScore& bls, const float& blsThreshold, const std::unordered_set<STPosition, STPositionHash>& positions,
    const std::string& currentMotif, const float blsScore,  int curDegenerateLetters, MotifCollection& processed, std::ostream& out) const
{
    // std::cerr << "currentMotif [" << currentMotif;
    // std::cerr << "currentMotif [" <<currentMotif << "\t" << blsScore << "]" << std::endl;
    if((unsigned char) currentMotif.length() >= l.first) {
        std::string representation = Motif::getRepresentative(currentMotif);
        if(processed.checkAndAddElement(representation)) {
            out << representation << "\t" << std::endl;
        }
    } // currentMotif is part of the printable
    // return if length is l.second - 1, since we cannot extend anymore then
    if((unsigned char) currentMotif.length() == l.second - 1) { return; }
    // std::cerr << " can be extended" << std::endl;
    // extend with possible characters in alphabet

    std::unordered_set<STPosition, STPositionHash> newpositions;
    int degenerateCount;
    std::bitset<N_BITS> occurence(0);
    const std::vector<IupacMask>* curalphabet = (curDegenerateLetters == maxDegenerateLetters) ?  &exactAlphabet : this->alphabet;
    for (IupacMask extension : *curalphabet) {
        // std::cerr << currentMotif << " + " << extension << std::endl;
        // increment position in positions list if possible
        degenerateCount = extension.isDegenerate() ? curDegenerateLetters + 1 : curDegenerateLetters;
        // if(degenerateCount > maxDegenerateLetters) { continue; } // cannot extend this degenerate character

        advanceIupacCharacter(extension, positions, newpositions, occurence);

        // can be extended if at least one new position is found!
        if(newpositions.size() > 0) {
            // std::cout << currentMotif << " -> " << extension << " has " << newpositions.size() << " paths : ";
            // for (auto tmpp : newpositions) { std::cout << " at " << tmpp.node->begin() << " + " << tmpp.offset << " - ";}
            // std::cout << std::endl;
            float childBlsScore = bls.getBLSScore(occurence);
            // std::cerr << "occ " << occurence << " for " << (currentMotif + extension.getRepresentation()) << " gives a bls score of " << childBlsScore << std::endl;
            if(childBlsScore >= blsThreshold) {
                recPrintMotifs(l, maxDegenerateLetters, bls, blsThreshold,
                               newpositions, currentMotif + extension.getRepresentation(), childBlsScore,
                               degenerateCount, processed, out);
            }
        }
    }
    newpositions.clear();
}

void SuffixTree::advanceIupacCharacter(IupacMask mask, const std::unordered_set<STPosition, STPositionHash>& currentPositions,
             std::unordered_set<STPosition, STPositionHash>& newPositions, std::bitset<N_BITS>& occurence) const {


    // std::cerr << "adding " << mask.getRepresentation() << std::endl;
    const std::vector<char> *chars = mask.getCharacters(); // chars is cleared in this function
    occurence.reset();
    newPositions.clear();
    for(STPosition p : currentPositions) {
        for(char c : *chars) { // these are the chars that belong to the iupac character extension!!
            if (p.atNode()) {
                // if p extended the full branch then go to child node with char c
                STNode *child = p.node->getChild(c);
                if(child != NULL) {
                    newPositions.insert(STPosition(child, 1));
                    occurence = occurence | child->getOccurence();
                    // std::cerr << "update occ " << occurence << std::endl;
                }
            } else {
                // case b) we are at an edge: try and match next character along edge
                if (T[p.node->begin() + p.offset] == c) {
                    // std::cerr << currentMotif << " extension to child " << c << " at " << p.node->begin() << " + " << p.offset<< std::endl;
                    occurence = occurence | p.node->getOccurence();
                    // std::cerr << "update occ " << occurence << std::endl;
                    newPositions.insert(STPosition(p.node, p.offset + 1));
                }
            }
        }
    }
}


// EXACT
const std::vector<IupacMask> SuffixTree::exactAlphabet ({
    IupacMask(BASE_A),
    IupacMask(BASE_C),
    IupacMask(BASE_G),
    IupacMask(BASE_T)
    });
// EXACT AND N
const std::vector<IupacMask> SuffixTree::exactAndNAlphabet ({
    IupacMask(BASE_A),
    IupacMask(BASE_C),
    IupacMask(BASE_G),
    IupacMask(BASE_T),
    IupacMask(IUPAC_N)
    });
// EXACT, TWOFOLD AND N
const std::vector<IupacMask> SuffixTree::exactTwofoldAndNAlphabet ({
    IupacMask(BASE_A),
    IupacMask(BASE_C),
    IupacMask(BASE_G),
    IupacMask(BASE_T),
    IupacMask(IUPAC_N),
    IupacMask(IUPAC_R),
    IupacMask(IUPAC_Y),
    IupacMask(IUPAC_S),
    IupacMask(IUPAC_W),
    IupacMask(IUPAC_K),
    IupacMask(IUPAC_M)
    });

// ============================================================================
// SUFFIX TREE (PUBLIC FUNCTIONS)
// ============================================================================


SuffixTree::SuffixTree(const string& T, bool hasReverseComplement) : T(T), reverseComplementFactor(hasReverseComplement ? 2 : 1)
{
        // maximum string length = 2^32-1
        if (T.size() >= (size_t)numeric_limits<length_t>::max())
                throw runtime_error("String exceeds maximum length");

        // construct suffix tree using Ukonen's algorithm
        constructUkonen();

        // construct suffix tree using naive algorithm
        //constructNaive();
}



SuffixTree::~SuffixTree()
{
        // Depth-first traversal of the tree
        stack<STNode*> stack;
        stack.push(root);

        while (!stack.empty()) {
                STNode* node = stack.top();
                stack.pop();

                for (int i = 0; i < MAX_CHAR; i++)
                        if (node->getChild(i) != NULL)
                                stack.push(node->getChild(i));

                delete node;
        }
}

void SuffixTree::printMotifs(const std::pair<short, short>& l, const Alphabet alphabet, const int& maxDegenerateLetters, const BLSScore& bls, const float& blsThreshold, std::ostream& out)
{
        if(alphabet == EXACT) {
            assert(maxDegenerateLetters == 0); // cannot have degenerate letters with exact alphabet
        }
        if(alphabet == EXACT)
            this->alphabet = &(this->exactAlphabet);
        else if(alphabet == EXACTANDN)
            this->alphabet = &(this->exactAndNAlphabet);
        else
            this->alphabet = &(this->exactTwofoldAndNAlphabet);

        MotifCollection processed;

        // std::string iupacword = "AAACMAKMTTT";
        // std::bitset<N_BITS> occurence(0);
        // std::unordered_set<STPosition, STPositionHash> positions;
        // matchIupacPattern(iupacword, positions, occurence);
        // if(!positions.empty()) {
        //     std::cerr << iupacword << " found with occurence " << occurence << " at " << positions.size() << " locations" << std::endl;
        //     for(auto p : positions) {
        //         std::cerr << "pos: " << p.getPositionInText() << ": " <<  T.substr(p.getPositionInText() - p.getDepth(), iupacword.length()) << std::endl;
        //     }
        // }
        // recPrintMotifs(l, maxDegenerateLetters, bls, blsThreshold, positions, iupacword, bls.getBLSScore(occurence), 3, processed, out);
        // std::cerr << "--------------------------" << std::endl;
        std::unordered_set<STPosition, STPositionHash> positions{STPosition(root)};
        recPrintMotifs(l, maxDegenerateLetters, bls, blsThreshold, positions, "", 1.0, 0, processed, out);
        std::cerr << "processed " << processed.size() << " motifs" << std::endl;
}


void SuffixTree::matchPattern(const string& P, vector<size_t>& occ)
{
        // clear the occurrence vector
        occ.clear();

        // can we completely match P?
        STPosition pos(root);
        if (!advancePos(pos, P, 0, P.size()))
                return;

        // get all occurrences under position
        getOccurrences(pos.node, occ);
}

void SuffixTree::matchIupacPattern(const string& P, std::unordered_set<STPosition, STPositionHash>& positions, std::bitset<N_BITS>& occurence)
{
        // can we completely match P?
        std::unordered_set<STPosition, STPositionHash> currentPositions{STPosition(root)};
        size_t i =0;

        while (!currentPositions.empty() && i < P.size()) {
            advanceIupacCharacter(IupacMask::characterToMask[P[i]], currentPositions, positions, occurence);
            currentPositions = positions;
            // std::cerr << "matching " << P[i] << " has " << currentPositions.size() << " locations and occurence " << occurence << std::endl;
            // for(auto p : positions) {
            //     std::cerr << "pos: " << p.getPositionInText() << ": " << T.substr(p.getPositionInText() - p.getDepth(), i + 1) << std::endl;
            // }
            i++;
        }
}


void SuffixTree::matchPattern(const string& P, BLSScore& bls)
{

        // can we completely match P?
        STPosition pos(root);
        if (!advancePos(pos, P, 0, P.size()))
                return;


        std::cerr << "Found occurrence @ " << pos.getPositionInText()
                  << " with occurence vector " << pos.getOccurence()
                  << " and BLS score " << bls.getBLSScore(pos.getOccurence())
                  << std::endl;
}

void SuffixTree::findMEM(const string& Q, size_t minSize, vector<MEMOcc>& occ)
{
        // clear the occurrence vector
        occ.clear();

        STPosition pos(root);
        for (size_t j = 0; j < Q.size(); j++) {

                // advance posMax as much as possible
                advancePos(pos, Q, j+pos.getDepth(), Q.size());

                if (pos.getDepth() >= minSize)
                        reportMEM(Q, j, minSize, pos, occ);

                pos = followSuffixLink(pos);
        }
}

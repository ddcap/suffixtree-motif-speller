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
#include <list>
#include "suffixtree.h"
#include "motif.h"

using namespace std;

// ============================================================================
// SUFFIX TREE (PRIVATE FUNCTIONS)
// ============================================================================

const std::vector<char> STNode::Alphabet ({ 'A', 'C', 'G', 'T', 'N', ' ', '$' });
// '$' is deliminter, ' ' is used to separate genes of same species and  'N'  is used to max exons
// const std::vector<short> STNode::charToIndex ({
const short STNode::charToIndex[] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 0
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 16
     5, -1, -1, -1,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 32
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 48
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1,  4, -1, // 64
    -1, -1, -1, -1,  3 // 80
};

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
        node_count++;
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
        node_count++;
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
        node_count++;


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
        node_count ++;

        // algorithm invariant: pos points to T[i:j-1[
        STPosition pos(root);
        unsigned char currentleafbit = 0; // supports up to 255 strings!
        int currentNumLeavesStart = 1;

        // in phase j, build implicit suffix tree for prefix T[0:j[
        for (size_t j = 1, numLeaves = 0; j <= T.size(); j++) {
                STNode *prevInternal = NULL;
                // i starts from numleaves -> keep the current bit to set in a variable!

                // skip 'numleafs' times rule 1 (extension of leaf)
                unsigned char currentbit = currentleafbit;
                int currentIStart = currentNumLeavesStart;
                for (size_t i = numLeaves; i < j; i++) {
                        // std::cout << "current i " << i << " is at " << T[i] << std::endl;
                        // note that pos will always point to T[i:j-1[ at this point
                        // if t[i] == # then set bit of currentbit + increment currentbit
                        // add bit to existing nodes


                        // add a SL from the previously created internal node
                        if (prevInternal != NULL && pos.atNode()) {
                                prevInternal->setSuffixLink(pos.node);
                                prevInternal = NULL;
                        }

                        // try and add character T[j-1]
                        if (advancePos(pos, T[j-1])) {
                            if(numLeaves != i && i + 1 == stringStartPositions[currentIStart]) { // only increment one if i == numleaves
                                // std::cout << "incrementing currentbit at " << T[i] << std::endl;
                                currentbit++;
                                currentIStart++;
                            } // increment currentbit for next iteration in for loop.
                            // std::cerr << "break at " << i << std::endl;
                            break; // rule 3 : do nothing + show stopper

                        }

                        // rule 2 : create internal node (optionally) and leaf
                        if (!pos.atNode()) {
                                pos = splitEdge(pos); // pos points to new node
                                if (prevInternal != NULL)
                                        prevInternal->setSuffixLink(pos.node);
                                prevInternal = pos.node;
                        }

                        addLeaf(pos, i, currentbit);
                        if(i + 1 == stringStartPositions[currentIStart]) {
                            // std::cout << "incrementing currentbit at " << T[i] << std::endl;
                            currentbit++;
                            currentIStart++;
                        } // increment currentbit for next iteration in for loop.

                        // extension is complete: follow suffix link
                        pos = followSuffixLink(pos);

                        numLeaves++;
                        if(numLeaves == stringStartPositions[currentNumLeavesStart]) {
                            // std::cout << "incrementing currentleafbit at " << T[numLeaves] << std::endl;
                            currentleafbit++;
                            currentNumLeavesStart++;
                        } // if previous character where the suffix starts is # then a new string has started so increment currentleafbit
                        // std::cerr << *this << std::endl;
                }
        }
        std::cerr << "ST of length "<< T.size() <<  ", memory usage: " <<  ((sizeof(SuffixTree) + sizeof(STNode) * node_count) / 1024 / 1024) << "MB" << std::endl;
}

// Routines to explore SuffixTree
void SuffixTree::printMotifString(const short& maxlen, const std::string& currentMotif, const BLSScore& bls, const occurence_bits& occurence, std::ostream& out) {
    // std::cerr << "nodes that match " << currentMotif << ":  with occ " << +occurence << " and blsScore: " << bls.getBLSScore(occurence) << std::endl;
    if(Motif::isRepresentative(currentMotif)) {
        Motif::writeMotif(currentMotif, out);
        // Motif::writeGroupIDAndMotif(currentMotif, out);
        out << "\t";
        bls.writeBLSVector(occurence, out);
        out << '\n'; // std::endl has a flushline which destroys performance.
        motifCount++;
    }
}
void SuffixTree::printMotifBinary(const short& maxlen, const std::string& currentMotif, const BLSScore& bls, const occurence_bits& occurence, std::ostream& out) {
    // std::cerr << "nodes that match " << currentMotif << ":  with occ " << +occurence << " and blsScore: " << bls.getBLSScore(occurence) << std::endl;
    if(Motif::isRepresentative(currentMotif)) {
        Motif::writeMotifInBinary(currentMotif, maxlen, out);
        // Motif::writeGroupIDAndMotifInBinary(currentMotif, out);
        bls.writeBLSVectorInBinary(occurence, out);
        motifCount++;
    }
}
/**
Recurse into the tree with degenerate letters, meaning we need a vector of positions we are currently at!
with 3 N's the most number of positions is 4*4*4= 64 positions to check.
To check validity of the motif with regard to the BLS score, we do the or operation on the mask of every position and then check this occurence bitset in the bls score function.
If score is equal to or more than the threshold then we keep the motif if its less we should stop the recursion!
*/
void SuffixTree::recPrintMotifs(const std::pair<short, short>& l,
    const int& maxDegenerateLetters, const BLSScore& bls, STPositionsPerLetter& matchingNodes,
    const std::string& prefix, int curDegenerateLetters, std::ostream& out)
{
    occurence_bits occurence(0);
    const std::vector<IupacMask>* curalphabet = (curDegenerateLetters == maxDegenerateLetters) ?  &exactAlphabet : this->alphabet;

    for (IupacMask extension : *curalphabet) {
        // increment position in positions list if possible
        if(extension.isDegenerate()){
            advanceIupacCharacter(extension, prefix.size(), matchingNodes, occurence);
        } else {
            advanceExactCharacter(extension, prefix.size(), matchingNodes, occurence);
        }

        // can be extended if at least one new position is found!
        if(matchingNodes.list[prefix.size() + 1].validPositions > 0) {
            std::string currentMotif = prefix + extension.getRepresentation();
            iteratorCount++;
            if(iteratorCount % 1000000 == 0) std::cerr << "\33[2K\r" << iteratorCount / 1000000 << " M motifs iterated" << std::flush;
            if(bls.greaterThanMinThreshold(occurence)) {

               if((unsigned char) currentMotif.length() >= l.first) { // print motif if correct length!
                    (this->*printMotif)(l.second, currentMotif, bls, occurence, out);
               }

                if((unsigned char) currentMotif.length() +  1 == l.second) { // max length reached, we do not recurse into the next extension
                    continue; // continue for loop/go to next extension
                } else { // recursive function to add one more letter to the motifs
                    recPrintMotifs(l, maxDegenerateLetters, bls, matchingNodes, currentMotif,
                                   (extension.isDegenerate() ? curDegenerateLetters + 1 : curDegenerateLetters), out);
               }
            }
        }
    }
}




// the return vector is a vector of positions, <# of string, pos in that string>
void SuffixTree::recPrintMotifsWithPositions(const std::pair<short, short>& l,
    const int& maxDegenerateLetters, const BLSScore& bls, STPositionsPerLetter& matchingNodes, std::vector<std::pair<int, int>>& stringPositions,
    const std::string& prefix, int curDegenerateLetters, std::ostream& out)
{
    occurence_bits occurence(0);

//TODO fix max positions -> but still find postions here and return to above
    int i = 0;
    int validChildren = 0;
    // std::string motif = "AAAAAARWCARA";

    if((unsigned char) prefix.length() +  1 < l.second) { // only extend if possible, else only get positions
        const std::vector<IupacMask>* curalphabet = (curDegenerateLetters == maxDegenerateLetters) ?  &exactAlphabet : this->alphabet;
        for (IupacMask extension : *curalphabet) {
            if(extension.isDegenerate()){
                advanceIupacCharacter(extension, prefix.size(), matchingNodes, occurence);
            } else {
                advanceExactCharacter(extension, prefix.size(), matchingNodes, occurence);
            }

            // if(motif == (prefix + extension.getRepresentation()).substr(0, motif.length())) {
            //     std::cerr << prefix + extension.getRepresentation() << " found at " << matchingNodes.list[prefix.size() + 1].validPositions << " positions\n";
            // }
            // can be extended if at least one new position is found!
            if(matchingNodes.list[prefix.size() + 1].validPositions > 0) {
                std::vector<std::pair<int, int>> newPositions;
                std::string currentMotif = prefix + extension.getRepresentation();
                if ((unsigned char) currentMotif.length() +  1 == l.second) {
                    getLeafPositions(newPositions, matchingNodes.list[prefix.size() + 1].list, matchingNodes.list[prefix.size() + 1].validPositions);
                } else {
                    recPrintMotifsWithPositions(l, maxDegenerateLetters, bls, matchingNodes, newPositions,
                                                currentMotif, (extension.isDegenerate() ? curDegenerateLetters + 1 : curDegenerateLetters), out);
                }
                if(!extension.isDegenerate()) { // only add positions of exact extensions
                    validChildren += matchingNodes.list[prefix.size() + 1].validPositions;
                    for (auto pos: newPositions) { // add to positions of this prefix
                         // std::cerr << "pos: " << pos.first << ", " << pos.second << std::endl;
                         stringPositions.push_back(pos);
                    }
                }

                if((unsigned char) currentMotif.length() >= l.first && newPositions.size() > 1) { // needs at least more than 1 read to have a bls score >0!!
                    getBestOccurence(newPositions, bls, occurence);
                    // if(motif == (currentMotif).substr(0, motif.length())) {
                    //     std::cerr << "valid positions found for " << currentMotif << " with " << newPositions.size() << " positions found with occ: " << +occurence << std::endl;
                    //     for (auto p: newPositions)
                    //         std::cerr << "pos : " << p.first <<", " << p.second << std::endl;
                    // }
                    if(bls.greaterThanMinThreshold(occurence)) { // print motif if correct length
                        (this->*printMotif)(l.second, currentMotif, bls, occurence, out);
                    }
                }
            }
            i++;
        }
    }
    if(validChildren == 0) {
        // std::cerr << prefix << " cannot be extended further, getting all positions from position list" << std::endl;
        getLeafPositions(stringPositions, matchingNodes.list[prefix.size()].list, matchingNodes.list[prefix.size()].validPositions);
    } else {
        // now find positions in current nodes that continue with a '-' (filler)
        getPositionsStartingWithDelimiter(stringPositions, matchingNodes.list[prefix.size()].list, matchingNodes.list[prefix.size()].validPositions);
    }

}

// TODO fix bug with scores here, is incorrect now...
void SuffixTree::getBestOccurence(std::vector<std::pair<int, int>>& positions, const BLSScore& bls, occurence_bits& occurence) {
    occurence = 0;
    float maxBls = 0;
    sort(positions.begin(), positions.end(), SuffixTree::positionPairSort);
    size_t i = 0;
    while(i < positions.size()) {
        occurence_bits motifOcc(0);
        occurence_bits rcOcc(0);
        int pos = positions[i].second;
        if((i + 1 < positions.size() && pos != positions[i].second) || i == positions.size() - 1) {
             // if only one occurence or new pos as last pos so also only one occurence increment and skip
            i++;
            continue;
        }
        while(i < positions.size() && pos == positions[i].second) {
            if(positions[i].first % reverseComplementFactor == 0) {
                motifOcc |= 1 << (positions[i].first / reverseComplementFactor);
            } else {
                rcOcc |= 1 << (positions[i].first / reverseComplementFactor);
            }
            i++;
        }
        // std::cerr << "pos: " << pos << " -> occ: " << +motifOcc << " rvOcc: " << +rcOcc << std::endl;
        if(__builtin_popcountll(motifOcc) > 1) {
            iteratorCount++;
            if(iteratorCount % 1000000 == 0) std::cerr << "\33[2K\r" << iteratorCount / 1000000 << " M motifs iterated" << std::flush;
            float score = bls.getBLSScore(motifOcc);
            if(score > maxBls) {
                maxBls = score;
                occurence = motifOcc;
            }
        }
        if(__builtin_popcountll(rcOcc) > 1) {
            iteratorCount++;
            if(iteratorCount % 1000000 == 0) std::cerr << "\33[2K\r" << iteratorCount / 1000000 << " M motifs iterated" << std::flush;
            float score = bls.getBLSScore(rcOcc);
            if(score > maxBls) {
                maxBls = score;
                occurence = rcOcc;
            }
        }
    }
}

void SuffixTree::getLeafPositions(std::vector<std::pair<int, int>>& positions, const std::vector<STPosition>& matchingNodes, const size_t size) const {
    std::vector<size_t> occ;
    for(size_t i = 0; i < size; i++) {
        getOccurrences(matchingNodes[i], occ);
    }
    std::sort(occ.begin(), occ.end());
    int stringId = 1;
    // std::cerr << occ.size() << " leaf positions from " << size << " node positions " << std::endl;
    for(auto p : occ) {
        while(p >= stringStartPositions[stringId]) { stringId++;} // find the correct string id for this occurence
        int posInString = p - stringStartPositions[stringId - 1];
        // std::cerr << p << " -> [" << stringId - 1 << ", " << posInString << "], ";
        positions.push_back(std::pair<int, int>(stringId - 1, posInString));
    }
    // std::cerr << std::endl;
}

void SuffixTree::getPositionsStartingWithDelimiter(std::vector<std::pair<int, int>>& positions, const std::vector<STPosition>& matchingNodes, const size_t size) const {
    std::vector<size_t> occ;
    for(size_t i = 0; i < size; i++) {
        // if node is not at
        if(matchingNodes[i].atNode()) {
            // if (matchingNodes[i].node->getChild(IupacMask::FILLER) != NULL)
                // getOccurrences(STPosition(matchingNodes[i].node->getChild(IupacMask::FILLER)), occ);
            // if (matchingNodes[i].node->getChild(IupacMask::STRINGDELIMITER) != NULL)
                // getOccurrences(STPosition(matchingNodes[i].node->getChild(IupacMask::STRINGDELIMITER)), occ);
            if (matchingNodes[i].node->getChild(IupacMask::DELIMITER) != NULL)
                getOccurrences(STPosition(matchingNodes[i].node->getChild(IupacMask::DELIMITER)), occ);
        } else {
             // if(T[matchingNodes[i].node->begin() + matchingNodes[i].offset] == IupacMask::FILLER ||
                // T[matchingNodes[i].node->begin() + matchingNodes[i].offset] == IupacMask::STRINGDELIMITER ||
            if(T[matchingNodes[i].node->begin() + matchingNodes[i].offset] == IupacMask::DELIMITER)
                getOccurrences(matchingNodes[i], occ);
        }
    }
    std::sort(occ.begin(), occ.end());
    int stringId = 1;
    // std::cerr << occ.size() << " '" << IupacMask::DELIMITER << "' positions from " << size << " node(s) ";
    for(auto p : occ) {
        while(p >= stringStartPositions[stringId]) { stringId++;} // find the correct string id for this occurence
        int posInString = p - stringStartPositions[stringId - 1];
        // std::cerr << p << " -> [" << stringId - 1 << ", " << posInString << "], ";
        positions.push_back(std::pair<int, int>(stringId - 1, posInString));
    }
    // std::cerr << std::endl;
}

void SuffixTree::advanceIupacCharacter(const IupacMask& mask, const int& characterPos, STPositionsPerLetter& positions, occurence_bits& occurence) const {
    // std::cerr << "adding " << mask.getRepresentation()  << " @ " << characterPos << std::endl;
    const std::string* chars = mask.getCharacters();
    occurence = 0; // reset!
    positions.list[characterPos + 1].reset();
    // std::cerr << "reset next letter list: " << positions.list[characterPos + 1].validPositions << std::endl;
    // std::cerr << "current positions list: " << positions.list[characterPos].validPositions << std::endl;
    for(size_t i = 0; i < positions.list[characterPos].validPositions; i++) {
    // for(STPosition p : currentPositions) {
        for(char c : *chars) { // these are the chars that belong to the iupac character extension!!
            if (positions.list[characterPos].list[i].atNode()) {
                // if p extended the full branch then go to child node with char c
                STNode *child = positions.list[characterPos].list[i].node->getChild(c);
                if(child != NULL) {
                    positions.list[characterPos + 1].addSTPosition(child, 1);
                    occurence |= child->getOccurence();
                    // std::cerr << "occ update " << +child->getOccurence() << " -> " << +occurence << std::endl;
                }
            } else {
                // case b) we are at an edge: try and match next character along edge
                if (T[positions.list[characterPos].list[i].node->begin() + positions.list[characterPos].list[i].offset] == c) {
                    positions.list[characterPos + 1].addSTPosition(positions.list[characterPos].list[i].node, positions.list[characterPos].list[i].offset + 1);
                    occurence |= positions.list[characterPos].list[i].node->getOccurence();
                    // std::cerr << "occ update " << +positions.list[characterPos].list[i].node->getOccurence() << " -> " << +occurence << std::endl;
                }
            }
        }
    }
    // std::cerr << "next letter list: " << positions.list[characterPos + 1].validPositions << std::endl;
}
void SuffixTree::advanceExactCharacter(const IupacMask& mask, const int& characterPos, STPositionsPerLetter& positions, occurence_bits& occurence) const {
    // std::cerr << "adding " << c  << " @ " << characterPos << std::endl;
    occurence = 0; // reset!
    const char c = mask.getRepresentation();
    positions.list[characterPos + 1].reset();
    for(size_t i = 0; i < positions.list[characterPos].validPositions; i++) {
        if (positions.list[characterPos].list[i].atNode()) {
            // if p extended the full branch then go to child node with char c
            STNode *child = positions.list[characterPos].list[i].node->getChild(c);
            if(child != NULL) {
                positions.list[characterPos + 1].addSTPosition(child, 1);
                occurence |= child->getOccurence();
                // std::cerr << "occ update " << +child->getOccurence() << " -> " << +occurence << std::endl;
            }
        } else {
            // case b) we are at an edge: try and match next character along edge
            if (T[positions.list[characterPos].list[i].node->begin() + positions.list[characterPos].list[i].offset] == c) {
                positions.list[characterPos + 1].addSTPosition(positions.list[characterPos].list[i].node, positions.list[characterPos].list[i].offset + 1);
                occurence |= positions.list[characterPos].list[i].node->getOccurence();
                // std::cerr << "occ update " << +positions.list[characterPos].list[i].node->getOccurence() << " -> " << +occurence << std::endl;
            }
        }
    }
    // std::cerr << "extension gives " << positions.list[characterPos + 1].validPositions << " valid positions" << std::endl;
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


SuffixTree::SuffixTree(const string& T, bool hasReverseComplement, std::vector<size_t> stringStartPositions_) :
    T(T), reverseComplementFactor(hasReverseComplement ? 2 : 1), stringStartPositions(stringStartPositions_)
{
        // for (size_t i = 0 ; i < stringStartPositions.size() - 1; i++) {
            // std::cerr << T.substr(stringStartPositions[i] , (stringStartPositions[i+1] - stringStartPositions[i])) << std::endl;
        // }
        // maximum string length = 2^32-1
        if (T.size() >= (size_t)numeric_limits<length_t>::max())
                throw runtime_error("String exceeds maximum length");

        // construct suffix tree using Ukonen's algorithm
        constructUkonen();

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

int SuffixTree::printMotifs(const std::pair<short, short>& l, const Alphabet alphabet, const int& maxDegenerateLetters, const BLSScore& bls, std::ostream& out, bool isAlignmentBased)
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

        STPositionsPerLetter positions(l.second, maxDegenerateLetters); // 13

// start from iupacword
        // std::string iupacword = "TNAGC";
        // occurence_bits occurence(0);
        // std::vector<STPosition> startPositions = matchIupacPattern(iupacword, 3, occurence);
        // for(auto p : startPositions) {
            // std::cerr << "pos: " << p.getPositionInText() << ": " <<  T.substr(p.getPositionInText() - p.getDepth(), iupacword.length()) << std::endl;
            // positions.list[iupacword.length()].addSTPosition(p.node, p.offset);
        // }
        // recPrintMotifs(l, maxDegenerateLetters, bls, positions, iupacword, 1, out);
// start from root
        motifCount = 0;
        iteratorCount = 0;
        positions.list[0].addSTPosition(root);
        std::vector<std::pair<int, int>> stringPos;
        if(isAlignmentBased)
            recPrintMotifsWithPositions(l, maxDegenerateLetters, bls, positions, stringPos, "", 0, out);
        else
            recPrintMotifs(l, maxDegenerateLetters, bls, positions, "", 0, out);
        return motifCount;
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

std::vector<std::pair<int, int>> SuffixTree::matchIupacPatternWithPositions(const string& P, const BLSScore& bls, int maxDegenerateLetters, occurence_bits& occurence)
{
        size_t i =0;

        STPositionsPerLetter positions(P.size() + 1, maxDegenerateLetters);
        positions.list[0].addSTPosition(root);
        while (!positions.list[i].empty() && i < P.size()) {
            // std::cerr << "valid pos: " << positions.list[i].validPositions << std::endl;
            advanceIupacCharacter(IupacMask::characterToMask[P[i]], i, positions, occurence);
            // std::cerr << "matching " << P[i] << " has " << positions.list[i+1].validPositions << " locations and occurence " << occurence << std::endl;
            // for(int j = 0; j < positions.list[i+1].validPositions; j++) {
            //     std::cerr << "pos: " << positions.list[i+1].list[j].getPositionInText() << ": " <<
            //     T.substr(positions.list[i+1].list[j].getPositionInText() - positions.list[i+1].list[j].getDepth(), i + 1) << std::endl;
            // }
            i++;
        }
        std::vector<std::pair<int, int>> stringPositions;
        // std::cerr << "found " << positions.list[P.size()].validPositions << " valid positions" << std::endl;
        // for(size_t i = 0; i < positions.list[P.size()].validPositions; i++) {
        //     std::cerr << i << ": " << positions.list[P.size()].list[i].getPositionInText() << std::endl;
        // }
        getLeafPositions(stringPositions, positions.list[P.size()].list, positions.list[P.size()].validPositions);
        // getPositionsStartingWithFiller(stringPositions, positions.list[P.size()].list, positions.list[P.size()].validPositions);
        getBestOccurence(stringPositions, bls, occurence);

        return stringPositions;
}


std::vector<STPosition> SuffixTree::matchIupacPattern(const string& P, int maxDegenerateLetters, occurence_bits& occurence)
{
        size_t i =0;

        STPositionsPerLetter positions(P.size() + 1, maxDegenerateLetters);
        positions.list[0].addSTPosition(root);
        while (!positions.list[i].empty() && i < P.size()) {
            // std::cerr << "valid pos: " << positions.list[i].validPositions << std::endl;
            advanceIupacCharacter(IupacMask::characterToMask[P[i]], i, positions, occurence);
            // std::cerr << "matching " << P[i] << " has " << positions.list[i+1].validPositions << " locations and occurence " << occurence << std::endl;
            // for(int j = 0; j < positions.list[i+1].validPositions; j++) {
                // std::cerr << "pos: " << positions.list[i+1].list[j].getPositionInText() << ": " <<
                // T.substr(positions.list[i+1].list[j].getPositionInText() - positions.list[i+1].list[j].getDepth(), i + 1) << std::endl;
            // }
            i++;
        }
        std::vector<STPosition> pos;
        for(size_t i = 0; i < positions.list[P.size()].validPositions; i++) {
            pos.push_back(positions.list[P.size()].list[i]);
        }
        return pos;
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



#include <vector>
#include <stack>
#include <iostream>
#include <bitset>
#include <random>
#include "motif.h"


//BLSLinkedListNode<N>
template<unsigned char N>
float BLSLinkedListNode<N>::getScore(const std::bitset<N>& occurence) {
    if(occurence.count() == 1) return 0.0f;
    // loop over next
    float score = 0.0f;
    int count = 0;
    std::vector<BLSLinkedListNode<N> *> list;
    BLSLinkedListNode<N> *iterator = this;
    while(iterator != NULL) {
        if((occurence & iterator->mask).any()) {
            count++;
            list.push_back(iterator);
            if(iterator->child == NULL) { // is a leaf add the score, independant of how many branches have occurence, since parent needs to be connected
                score += iterator->length;
            }
        }
        iterator = iterator->next;
    }
    if(count == 1) {
        if(list.front()->child != NULL) { // has children
            score += list.front()->child->getScore(occurence); // add the children since this branch has more than 1 edges that need to be connected
        }
    } else if(count > 1) {
        for(auto node : list) {
            if(node->child != NULL) { // has children
                score += node->length; // add length only if it hasnt been added before
                score += node->child->getScore(occurence); // add the children since this branch has more than 1 edges that need to be connected
            }
        }
    }

    return score;
}

template<unsigned char N>
std::ostream& BLSLinkedListNode<N>::write(std::ostream& o) const {
    for(int i = 0; i < level; i++ ) { o << "  "; }
    o << "[" << mask << "/" << length << "]"; // print this actual node
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


// BLSSCORE
void BLSScore::recReadBranch(int recursion, int& leafcount, std::string& newick, BLSLinkedListNode<N_BITS>* currentroot) {

    // int startleafcount = leafcount;
    BLSLinkedListNode<N_BITS>* currentnode = currentroot;
    bool endofbranch = false;
    // std::cout << "rec[" << recursion << "] starts at " << leafcount << std::endl;
    // std::cout << "rec[" << recursion << "] node: " << currentnode << std::endl;
    while(!newick.empty() && !endofbranch) {
        if(newick[0] == ';') { // end of line!
            // std::cout << "rec[" << recursion << "] ends all!" << std::endl;
            newick.erase(0,1);
        } else if(newick[0] == '(') { // new branch
            // add new branch start;
            newick.erase(0,1);
            BLSLinkedListNode<N_BITS>* child = currentnode->addChild(recursion + 1);
            recReadBranch(recursion + 1, leafcount, newick, child);
            // std::cout << "rec[" << recursion << "] received children " << std::endl << *child << std::endl;
            BLSLinkedListNode<N_BITS>* iterator = child;
            std::bitset<N_BITS> mask;
            while(iterator != NULL) {
                // std::cout << "rec[" << recursion << "] child mask: " << iterator->getMask() << std::endl;
                mask |= iterator->getMask();
                iterator = iterator->getNext();
            }
            currentnode->setMask(mask);
            // std::cout << "rec[" << recursion << "] read branch with mask " << mask << std::endl;
        } else if (newick[0] == ')') {
            newick.erase(0,1);
            endofbranch = true;
        } else if (newick[0] == ','){
            currentnode = currentnode->addNext(recursion);
            newick.erase(0,1); // just go to next leaf/branch
        } else if(newick[0] == ':') { // get length of branch and either add a node or a length to the list of branches
            newick.erase(0,1);
            // read score
            int charPosAfterScore = newick.find_first_not_of(".0123456789E-"); // E- for scientific score
            float length = std::stof(newick.substr(0, charPosAfterScore));
            newick.erase(0, charPosAfterScore);
            currentnode->setLength(length);
            // std::cout << "rec[" << recursion << "] read length: " << length << std::endl;
        } else {
            // read name + score
            int doublepointposition = newick.find_first_of(':');
            std::string name = newick.substr(0, doublepointposition);
            newick.erase(0, doublepointposition);
            std::bitset<N_BITS> mask;
            mask.set(leafcount);
            // std::cout << "rec[" << recursion << "] reading leaf " << name << " @ " << leafcount  << " " << mask << std::endl;
            currentnode->setMask(mask);
            leafcount++;

        }
    }
    // std::cout << "rec[" << recursion << "] ends at " << leafcount << std::endl;
}


float BLSScore::getBLSScore(std::bitset<N_BITS> occurence) const {
    return root->getChild()->getScore(occurence); // root has no next but only children! so first branch is of length 0 with 11111 so always true!
}

// IUPACMASK
void IupacMask::getCharacters(std::vector<char>& chars) {
    chars.clear();
    if (0x1 & mask) { chars.push_back('A'); }
    if (0x2 & mask) { chars.push_back('C'); }
    if (0x4 & mask) { chars.push_back('G'); }
    if (0x8 & mask) { chars.push_back('T'); }
}

char IupacMask::getRepresentation() const {
    switch (mask) {
    case BASE_A:  return 'A';
    case BASE_C:  return 'C';
    case IUPAC_M: return 'M';
    case BASE_G:  return 'G';
    case IUPAC_R: return 'R';
    case IUPAC_S: return 'S';
    case IUPAC_V: return 'V';
    case BASE_T:  return 'T';
    case IUPAC_W: return 'W';
    case IUPAC_Y: return 'Y';
    case IUPAC_H: return 'H';
    case IUPAC_K: return 'K';
    case IUPAC_D: return 'D';
    case IUPAC_B: return 'B';
    case IUPAC_N: return 'N';
    default: return 'N';
    }
}

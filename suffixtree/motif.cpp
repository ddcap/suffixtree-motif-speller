

#include <vector>
#include <stack>
#include <iostream>
#include <bitset>
#include <random>
#include "motif.h"

//MOTIF

const std::vector<unsigned char> Motif::complement ({
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 16
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 32
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 48
    0, 'T', 'V', 'G', 'H', 0, 0, 'C', 'D', 0, 0, 'M', 0, 'K', 'N', 0, // 64
    0, 0, 'Y', 'S', 'A', 0, 'B', 'W', 0, 'R' // 80
});

std::string Motif::ReverseComplement(std::string read) {
    std::string rc = "";
    for(int i = read.length() - 1; i>=0; i--) {
        rc += complement[read[i]];
    }
    return rc;
}

bool Motif::isRepresentative(std::string read) {
    size_t i = 0;
    while (i < read.length() && read[i] == complement[read[read.length() - 1 - i]]) {
        i++;
    }
    return read[i] < complement[read[read.length() - 1 - i]];
}

std::string Motif::getRepresentative(std::string read) {
    size_t i = 0;
    std::string rc = "";
    char rci = complement[read[read.length() - 1 - i]];
    while (i < read.length() && read[i] == rci) {
        rc += rci;
        i++;
        rci = complement[read[read.length() - 1 - i]];
    }
    if (read[i] < rci) {
        return read;
    } else {
        // return rc!!
        while (i < read.length()) {
            rc += rci;
            i++;
            rci = complement[read[read.length() - 1 - i]];
        }
        return rc;
    }
}

//BLSLinkedListNode<N>
template<unsigned char N>
float BLSLinkedListNode<N>::getScore(const std::bitset<N>& occurence) {
    if(occurence.count() <= 1) return 0.0f;
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

// template<unsigned char N>
// std::ostream& BLSLinkedListNode<N>::write(std::ostream& o) const {
//     for(int i = 0; i < level; i++ ) { o << "  "; }
//     o << "[" << mask << "/" << length << "]"; // print this actual node
//     if(child != NULL) {
//         o << std::endl;
//         o << *child;
//     }
//     if(next != NULL) {
//         o << std::endl;
//         o << *next;
//     }
//     return o;
// }


// BLSSCORE
void BLSScore::prepAllCombinations() {
    // precalculate all combinations of occurence
    for (int i = 0; i < pow(2, N_BITS); i++) {
        std::bitset<N_BITS> occurence(i);
        preparedBLS.push_back(calculateBLSScore(occurence));
        // std::cerr << "prepping " << occurence << " => " << occurence.to_ulong() << " " << preparedBLS[i] << std::endl;
    }
}

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


float BLSScore::calculateBLSScore(std::bitset<N_BITS> occurence) const {
    return root->getChild()->getScore(occurence); // root has no next but only children! so first branch is of length 0 with 11111 so always true!
}
float BLSScore::getBLSScore(std::bitset<N_BITS> occurence) const {
    return preparedBLS[occurence.to_ulong()];
}

// IUPACMASK
const std::vector<char>* IupacMask::getCharacters() {
    return &characterLists[mask];
}
const std::vector<IupacMask> IupacMask::characterToMask ({
    IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), // 0
    IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), // 16
    IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), // 32
    IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), IupacMask(), // 48
    IupacMask(), IupacMask(BASE_A), IupacMask(IUPAC_B), IupacMask(BASE_C), IupacMask(IUPAC_D), IupacMask(), IupacMask(), IupacMask(BASE_G), IupacMask(IUPAC_H), IupacMask(), IupacMask(), IupacMask(IUPAC_K), IupacMask(), IupacMask(IUPAC_M), IupacMask(IUPAC_N), IupacMask(), // 64
    IupacMask(), IupacMask(), IupacMask(IUPAC_R), IupacMask(IUPAC_S), IupacMask(BASE_T), IupacMask(), IupacMask(IUPAC_V), IupacMask(IUPAC_W), IupacMask(), IupacMask(IUPAC_Y) // 80
});
const std::vector<std::vector<char>> IupacMask::characterLists (
{
    {}, // 0!!
    {'A'},
    {'C'},
    {'A', 'C'},
    {'G'},
    {'A', 'G'},
    {'G', 'C'},
    {'A', 'C', 'G'},
    {'T'},
    {'A', 'T'},
    {'C', 'T'},
    {'A', 'C', 'T'},
    {'G', 'T'},
    {'A', 'G', 'T'},
    {'C', 'G', 'T'},
    {'A', 'C', 'G', 'T'},
});
const std::vector<char> IupacMask::representation ({ 0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N' });

char IupacMask::getRepresentation() const {
    return representation[mask];
}

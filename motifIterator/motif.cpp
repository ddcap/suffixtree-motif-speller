

#include "motif.h"
// #include <bitset>

//MOTIF

const std::vector<char> Motif::complement ({
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 16
    ' ', 0, 0, '#', '$', 0, 0, 0, 0, 0, 0, 0, 0, '-', 0, 0, // 32 // this is for the alignment based RC, replace '-' with '-'
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 48
    0, 'T', 'V', 'G', 'H', 0, 0, 'C', 'D', 0, 0, 'M', 0, 'K', 'N', 0, // 64
    0, 0, 'Y', 'S', 'A', 0, 'B', 'W', 0, 'R' // 80
});

std::string Motif::ReverseComplement(const std::string& read) {
    std::string rc = "";
    for(int i = read.length() - 1; i>=0; i--) {
        rc += complement[read[i]];
    }
    return rc;
}

std::string Motif::getGroupID(const std::string& read) {
      std::string tmp = read;
      std::sort(tmp.begin(), tmp.end());
      return tmp;
}

bool Motif::isRepresentative(const std::string& read) {
    size_t i = 0;
    while (i < read.length() && read[i] == complement[read[read.length() - 1 - i]]) {
        i++;
    }
    // if the same, ie i == read.length -> RC wont be matched in the tree since it is the same!!1 so no need to save them to check if it already has passed in the motifs.
    return i== read.length() || read[i] < complement[read[read.length() - 1 - i]];
}
bool Motif::isGroupRepresentative(const std::string& read) {
    std::string group = getGroupID(read);
    size_t i = 0;
    while (i < group.length() && group[i] == complement[group[group.length() - 1 - i]]) {
        i++;
    }
    // if the same, ie i == read.length -> RC wont be matched in the tree since it is the same!!1 so no need to save them to check if it already has passed in the motifs.
    return i== group.length() || group[i] < complement[group[group.length() - 1 - i]];
}

void Motif::writeMotif(const std::string& motif, std::ostream& out) {
    out << motif;
}
void Motif::writeGroupIDAndMotif(const std::string& motif, std::ostream& out) {
    out << getGroupID(motif) << "\t" << motif;
}
void Motif::writeMotifInBinary(const std::string& motif, const short &maxlen, std::ostream& out) {
    char size = motif.length(); // assumes length isnt more than 255 chars
    // write size
    out.write(&size, 1);
    char numberOfBytes = maxlen >> 1; // works since this maxlen is non inclusive (< instead of <=)

    //write Motif
    for(int i = 0; i < numberOfBytes; i++) {
        char toWrite = 0;
        if(i*2 < size)
            toWrite |= IupacMask::characterToMask[motif[i*2]].getMask() << 4;
        if(i*2 + 1 < size)
            toWrite |= IupacMask::characterToMask[motif[i*2+1]].getMask(); // << 4; // makes more sense if we combine it into a long later
        out.write(&toWrite, 1);
    }
}

void Motif::writeGroupIDAndMotifInBinary(const std::string& motif, const short &maxlen, std::ostream& out) {
    char size = motif.length(); // assumes length isnt more than 255 chars
    // write size
    out.write(&size, 1);
    char numberOfBytes = maxlen >> 1; // works since this maxlen is non inclusive (< instead of <=)

    std::string groupId = getGroupID(motif);
    for(int i = 0; i < numberOfBytes; i++) {
        char toWrite = 0;
        if(i*2 < size)
            toWrite |= IupacMask::characterToMask[groupId[i*2]].getMask() << 4;
        if(i*2 + 1 < size)
            toWrite |= IupacMask::characterToMask[groupId[i*2+1]].getMask(); // << 4;
        out.write(&toWrite, 1);
    }
    //write Motif
    for(int i = 0; i < numberOfBytes; i++) {
        char toWrite = 0;
        if(i*2 < size)
            toWrite |= IupacMask::characterToMask[motif[i*2]].getMask() << 4;
        if(i*2 + 1 < size)
            toWrite |= IupacMask::characterToMask[motif[i*2+1]].getMask(); // << 4;
        out.write(&toWrite, 1);
    }
}

std::string Motif::getRepresentative(const std::string& read) {
    size_t i = 0;
    std::string rc = "";
    char rci = complement[read[read.length() - 1 - i]];
    while (i < read.length() && read[i] == rci) {
        rc += rci;
        i++;
        rci = complement[read[read.length() - 1 - i]];
    }
    if (i < read.length() && read[i] < rci) {
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

//BLSLinkedListNode
float BLSLinkedListNode::getScore(const occurence_bits& occurence) {
    if(__builtin_popcountll(occurence) <= 1) return 0.0f;
    // loop over next
    float score = 0.0f;
    int count = 0;
    std::vector<BLSLinkedListNode *> list;
    BLSLinkedListNode *iterator = this;
    while(iterator != NULL) {
        if((occurence & iterator->mask) > 0) {
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
            float childscore = list.front()->child->getScore(occurence);
            score += childscore; // add the children since this branch has more than 1 edges that need to be connected
        }
    } else if(count > 1) {
        for(auto node : list) {
            if(node->child != NULL) { // has children
                score += node->length; // add length only if it hasnt been added before
                float childscore = node->child->getScore(occurence);
                score += childscore; // add the children since this branch has more than 1 edges that need to be connected
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
void BLSScore::prepAllCombinations(int used_bits) {
    // precalculate all combinations of occurence
    for (int i = 0; i < pow(2, used_bits); i++) {
        occurence_bits occurence(i);
        preparedBLS.push_back(calculateBLSScore(occurence));
        // std::cerr << std::bitset<16>(occurence) << "\t" << calculateBLSScore(occurence) << "\n";
        preparedBLSVector.push_back(calculateBLSVector(preparedBLS[i]));
    }
}

void BLSScore::recReadBranch(int recursion, int& leafcount, std::string& newick, BLSLinkedListNode* currentroot, std::vector<std::string> &order_of_species) {

    // int startleafcount = leafcount;
    BLSLinkedListNode* currentnode = currentroot;
    bool endofbranch = false;
    // std::cerr << "rec[" << recursion << "] starts at " << leafcount << std::endl;
    // std::cerr << "rec[" << recursion << "] node: " << currentnode << std::endl;
    while(!newick.empty() && !endofbranch) {
        if(newick[0] == ';') { // end of line!
            // std::cout << "rec[" << recursion << "] ends all!" << std::endl;
            newick.erase(0,1);
        } else if(newick[0] == '(') { // new branch
            // add new branch start;
            newick.erase(0,1);
            BLSLinkedListNode* child = currentnode->addChild(recursion + 1);
            recReadBranch(recursion + 1, leafcount, newick, child, order_of_species);
            BLSLinkedListNode* iterator = child;
            occurence_bits mask(0);
            int childcount = 0;
            while(iterator != NULL) {
                childcount++;
                // std::cerr << "rec[" << recursion << "] child mask: " << +iterator->getMask() << std::endl;
                mask |= iterator->getMask();
                iterator = iterator->getNext();
            }
            if(child->getNext() == NULL) { // merging single child
                currentnode->setChild(child->getChild());
                currentnode->addLength(child->getLength());
                delete child;
            }
            currentnode->setMask(mask);
            // std::cerr << "rec[" << recursion << "] read branch with mask " << +mask << std::endl;
        } else if (newick[0] == ')') {
            newick.erase(0,1);
            endofbranch = true;
        } else if (newick[0] == ','){
            currentnode = currentnode->addNext(recursion);
            newick.erase(0,1); // just go to next leaf/branch
        } else if(newick[0] == ':') { // get length of branch and either add a node or a length to the list of branches
            newick.erase(0,1);
            // read score
            int charPosAfterScore = newick.find_first_not_of(".0123456789Ee-"); // E- for scientific score
            float length = std::stof(newick.substr(0, charPosAfterScore));
            newick.erase(0, charPosAfterScore);
            currentnode->addLength(length);
            // std::cerr << "rec[" << recursion << "] read length: " << length << std::endl;
        } else {
            // read name + score
            int doublepointposition = newick.find_first_of(':');
            std::string name = newick.substr(0, doublepointposition);
            newick.erase(0, doublepointposition);
            occurence_bits mask(0);
            mask |= 1 << leafcount;
            order_of_species.push_back(name);
            // std::cout << "rec[" << recursion << "] reading leaf " << name << " @ " << leafcount  << " " << +mask << std::endl;
            currentnode->setMask(mask);
            leafcount++;

        }
    }
    // std::cout << "rec[" << recursion << "] ends at " << leafcount << std::endl;
}


// const std::vector<float> BLSScore::blsThresholds ({ 0.15, 0.5, 0.6, 0.7, 0.9, 0.95});

float BLSScore::calculateBLSScore(const occurence_bits& occurence) const {
    return root->getChild()->getScore(occurence); // root has no next but only children! so first branch is of length 0 with 11111 so always true!
}
float BLSScore::getBLSScore(const occurence_bits& occurence) const {
    return preparedBLS[occurence];
}
// TODO dont use vectors anymore but use unsigned char, indicating how many 1's -> how many BLS thresholds are met, assuming they are in ascending order!
char BLSScore::calculateBLSVector(const float& bls) const {
    std::vector<int> ret;
    char i = 1;
    while(i < blsThresholds.size() && bls > blsThresholds[i]) {
        i++;
    }
    return i;
}

// std::vector<int> BLSScore::calculateBLSVector(const float& bls) const {
//     std::vector<int> ret;
//     ret.push_back(1); // should have already past first check!
//     size_t i = 1;
//     while(i < blsThresholds.size() && bls > blsThresholds[i]) {
//         ret.push_back(1);
//         i++;
//     }
//     for(; i < blsThresholds.size(); i++) {
//         ret.push_back(0);
//     }
//     return ret;
// }
// const std::vector<int>* BLSScore::getBLSVector(const occurence_bits& occurence) const {
const char* BLSScore::getBLSVector(const occurence_bits& occurence) const {
    return &preparedBLSVector[occurence];
}

void BLSScore::writeBLSVector(const occurence_bits& occurence, std::ostream& out) const {
    //write Vector
    out << preparedBLSVector[occurence]; // [0];
    // for (size_t i = 1; i < preparedBLSVector[occurence].size(); i++) {
        // out << "," << preparedBLSVector[occurence][i];
    // }
}
void BLSScore::writeBLSVectorInBinary(const occurence_bits& occurence, std::ostream& out) const {
    // assume this is only 1 byte (max 8 thresholds), and first bit is first threshold and so on...
    // char size = preparedBLSVector[occurence].size();
    // out.write(&size, 1);

// TEST : print unsigned char of number of thresholds that are 1 -> one byte (unsigned) can have up to 256 thresholds!!!
    out.write(&preparedBLSVector[occurence], 1);


    //write Vector
    // char toWrite = 0; // convert to unsigned short int for more thresholds -> 16 thresholds
    // int i = 0;
    // for(; i < size; i++) {
    //     toWrite |= preparedBLSVector[occurence][i] << i;
    //     // if(i%8 == 7) { // write byte and reset for next
    //         // out.write(&toWrite, 1);
    //         // toWrite = 0;
    //     // }
    // }
    // // if(i%8 != 7)
    // out.write(&toWrite, 1);
}
char BLSScore::readBLSVectorInBinary(std::istream& in) const {
    char readChar = 0;
    in.read(&readChar, 1);
    return readChar;
}

bool BLSScore::greaterThanMinThreshold(const occurence_bits& occurence) const {
    return preparedBLS[occurence] > blsThresholds[0];
}

bool BLSScore::greaterThanThreshold(const occurence_bits& occurence, const int& blsThresholdIdx) const {
    return blsThresholdIdx < blsThresholds.size() && preparedBLS[occurence] > blsThresholds[blsThresholdIdx];
}

// IUPACMASK
char IupacMask::getRandomChar(char c) {
    const std::string *iupac_character_list = characterToMask[c].getCharacters();
    char new_c = iupac_character_list->at(rand() % iupac_character_list->size());
    std::cerr << "replaced " << c << " with " << new_c << std::endl;
    return new_c;
}

const std::string* IupacMask::getCharacters() const {
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
const std::vector<std::string> IupacMask::characterLists (
{
    {}, // 0!!
    {"A"},
    {"C"},
    {"AC"},
    {"G"},
    {"AG"},
    {"GC"},
    {"ACG"},
    {"T"},
    {"AT"},
    {"CT"},
    {"ACT"},
    {"GT"},
    {"AGT"},
    {"CGT"},
    {"ACGT"},
});
const std::vector<char> IupacMask::representation ({ 0, 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N' });

char IupacMask::getRepresentation() const {
    return representation[mask];
}

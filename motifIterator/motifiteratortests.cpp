#include "gtest/gtest.h"
#include <bits/stdc++.h>
#include <vector>
#include <algorithm>
#include "motif.h"
#include "suffixtree.h"


class MotifIteratorTest: public ::testing::Test {
protected:
    std::vector<float> blsThresholds{0.15, 0.5, 0.6, 0.7, 0.9, 0.95};
    std::pair<short, short> l;
    occurence_bits species_BD = 1;
    occurence_bits species_OS = 2;
    occurence_bits species_SB = 4;
    occurence_bits species_ZM = 8;
    int N;
    std::string name;
    BLSScore *bls = NULL;
    std::vector<std::string> order_of_species;
    std::vector<std::string> final_order_of_species{"BD", "OS", "SB", "ZM"};
    std::string newick = "((BD:0.2688,OS:0.2688):0.0538,(SB:0.086,ZM:0.086):0.2366);";
    std::string BD = "ACGACGTACGTACG"; // one matching motif of length 8: ACGTACGT
    std::string OS = "GCTACGTACGTGCT";
    std::string SB = "AGTACGTACGTAGT";
    std::string ZM = "TCTACGTACGTTCT";
    std::string BDAndRC = "ACGACGTACGTACG$CGTACGTACGTCGT$";
    std::string OSAndRC = "GCTACGTACGTGCT$AGCACGTACGTAGC$";
    std::string SBAndRC = "AGTACGTACGTAGT$ACTACGTACGTACT$";
    std::string ZMAndRC = "TCTACGTACGTTCT$AGAACGTACGTAGA$";
    std::vector<std::string> possible_motifs{"AACGTACG","ACGACGTA", "ACTACGTA", "AGAACGTA", "AGCACGTA", "AGTACGTA", "CACGTACG", "CGACGTAC", "CGTACGTC", "CGTACGTG", "CGTACGTT", "GAACGTAC", "GACGTACG", "GCACGTAC", "GCTACGTA", "GTACGTCG", "GTACGTGC", "GTACGTTC", "TACGTACT", "TACGTAGA", "TACGTAGC", "TACGTAGT", "TACGTCGT", "TACGTGCT", "TACGTTCT", "TCTACGTA"};
    std::vector<std::string> orthogroup{"TESTORTHO1", newick, "4", "BD\tBD", BD, "OS\tOS", OS, "ZM\tZM", ZM , "SB\tSB", SB};
    std::string finalT = "ACGACGTACGTACG$CGTACGTACGTCGT$GCTACGTACGTGCT$AGCACGTACGTAGC$TCTACGTACGTTCT$AGAACGTACGTAGA$AGTACGTACGTAGT$ACTACGTACGTACT$";
    std::string T;
    SuffixTree *ST = NULL;
    SuffixTree *STCounted = NULL;
    std::vector<size_t> stringStartPositions;
    std::vector<size_t> finalStringStartPositions{0, 15, 30, 45, 60, 75, 90, 105, 120};
    std::vector<size_t> order_of_species_mapping;
    std::vector<size_t> final_order_of_species_mapping{0, 1, 3, 2};
    std::vector<size_t> next_gene_locations;// only needed for find motifs, is same as startpositions of only 1 gene per species (no paralogues)
    std::vector<std::string> gene_names;
    size_t totalCount = 0;
    MyMotifMap *motif_to_blsvector_map = NULL;

    MotifIteratorTest( ) {
    }

    size_t getIndexOfVector(const std::vector<std::string> &v, const std::string &val) {
        auto it = find(v.begin(), v.end(), val);

        // If element was found
        int index = -1;
        if (it != v.end())
        {
            index = it - v.begin();
        } else {
            std::cerr << "Species not found in bls tree, will crash, please provide the correct names in the tree and the sequences below. " << val << std::endl;
        }
        return index;
    }

    void readOrthoGroup() {

        char blsvectorsize = (unsigned char)blsThresholds.size(); // assume its less than 256
        int i = 0;
        stringStartPositions.push_back(0);
        // READ DATA
        std::string line;
        name = orthogroup[i];
        i++;
        newick = orthogroup[i];
        i++;
        N = std::stoi(orthogroup[i]);
        bls = new BLSScore(blsThresholds, newick, N, order_of_species);
        i++;
        size_t current_pos = 0;
        next_gene_locations.push_back(current_pos);
        for (int j = 0; j < N; j++) {
            line = orthogroup[i];
            // gene names
            std::vector<std::string> genes;
            std::string species = line.substr(line.find_first_of('\t')+1);
            // std::cerr << species << std::endl;
            order_of_species_mapping.push_back(getIndexOfVector(order_of_species, species));
            line = line.substr(0, line.find_first_of('\t'));
            size_t start = 0;
            size_t end = line.find_first_of(' ', start);
            while(end != std::string::npos) {
                genes.push_back(line.substr(start, end - start));
                start = end + 1;
                end = line.find_first_of(' ', start);
            }
            genes.push_back(line.substr(start));
            for (size_t k =0; k < genes.size(); k++) {
                gene_names.push_back(genes[k]);
            } // add RC genes
            for (size_t k = genes.size() ; k > 0; k--) {
                gene_names.push_back(genes[k-1]);
            }
            // genes
            i++;
            line = orthogroup[i];
            if (!T.empty())
                T.push_back(IupacMask::DELIMITER);
            std::for_each(line.begin(), line.end(), [](char & c) { // convert all to upper case!
                if(c == IupacMask::FILLER)
                    c = IupacMask::DELIMITER;
                else {
                    c = ::toupper(c);
                    // should already be fixed in the preprocessing tool!!!!
                    // if (validCharacters.find(c) == validCharacters.end()) {
                        // c = IupacMask::getRandomChar(c);
                    // }
                }
            });
            T.append(line);
            T.push_back(IupacMask::DELIMITER);
            stringStartPositions.push_back(T.size());
            T.append(Motif::ReverseComplement(line));
            stringStartPositions.push_back(T.size() + 1);
            // std::cout << T << std::endl;

            // add gene start locations...
            std::vector<size_t> gene_sizes;
            start = 0;
            end = line.find_first_of(' ', start);
            while(end != std::string::npos) {
                gene_sizes.push_back(end + 1 - start);
                start = end + 1;
                end = line.find_first_of(' ', start);
            }
            gene_sizes.push_back(line.size() + 1 - start);
            for (size_t k =0; k < gene_sizes.size(); k++) {
                current_pos += gene_sizes[k];
                next_gene_locations.push_back(current_pos);
            } // add RC genes
            for (size_t k = gene_sizes.size(); k > 0; k--) {
                current_pos += gene_sizes[k - 1];
                next_gene_locations.push_back(current_pos);
            }
            i++;
        }
        T.push_back(IupacMask::DELIMITER);

        motif_to_blsvector_map = new MyMotifMap(blsvectorsize, l);
        ST = new SuffixTree(T, name, true, stringStartPositions, gene_names, next_gene_locations, order_of_species_mapping, NULL);
        STCounted = new SuffixTree(T, name, true, stringStartPositions, gene_names, next_gene_locations, order_of_species_mapping, motif_to_blsvector_map);

    }

    void SetUp( ) {
        // code here will execute just before the test ensues
        l = std::pair<short, short>(8, 9);
        readOrthoGroup();
    }

    void TearDown( ) {
        // code here will be called just after the test completes
        // ok to through exceptions from here if need be
        delete bls;
        delete ST;
        delete STCounted;
        delete motif_to_blsvector_map; // dont do this when i recursively print and delete!!
    }

    ~MotifIteratorTest( )  {
        // cleanup any pending stuff, but no exceptions allowed
    }
    // put in any custom data members that you need
};

std::map<char, char> complement {
    {'A', 'T'},
    {'C', 'G'},
    {'G', 'C'},
    {'T', 'A'}
};
std::string twofoldandN = "ACGTRYSWKMN";
std::map<char, std::string> expandchar{
    {'A', "A"},
    {'C', "C"},
    {'G', "G"},
    {'T', "T"},
    {'N', "ACGT"},
    {'R', "AG"},
    {'Y', "CT"},
    {'S', "GC"},
    {'W', "AT"},
    {'K', "GT"},
    {'M', "AC"},
    {'B', "CGT"},
    {'D', "AGT"},
    {'H', "ACT"},
    {'V', "ACG"},
};
std::string reversecomplement(std::string motif) {
    std::string rc = "";
    for (int i = 0; i < motif.length(); i++) {
        rc += complement[motif[motif.length() - 1 - i]];
    }
    return rc;
}
std::vector<std::string> expand(std::string motif) {
    std::vector<std::string> currentlist{""};
    for (int i = 0; i < motif.length(); i++) {
        std::vector<std::string> tmp;
        std::string expand = expandchar[motif[i]];
        for (int j = 0; j < currentlist.size(); j ++) {
            for (int k = 0; k < expand.length(); k ++) {
                tmp.push_back(currentlist[j] + expand[k]);
            }
        }
        currentlist = tmp;
    }
    // Add RC
    int listlen = currentlist.size();
    for (int i = 0; i < listlen; i++) {
        currentlist.push_back(reversecomplement(currentlist[i]));
    }
    return currentlist;
}

TEST_F (MotifIteratorTest, Setup) {
    ASSERT_EQ (order_of_species.size() , 4);
    ASSERT_EQ(final_order_of_species.size(), order_of_species.size());
    for(int i = 0; i < order_of_species.size(); i++ ) {
        ASSERT_EQ(final_order_of_species[i], order_of_species[i]);
    }
    ASSERT_EQ(T, finalT);
    ASSERT_EQ(finalStringStartPositions.size(), stringStartPositions.size());
    for(int i = 0; i < stringStartPositions.size(); i++ ) {
        ASSERT_EQ(finalStringStartPositions[i], stringStartPositions[i]);
    }
    ASSERT_EQ(final_order_of_species_mapping.size(), order_of_species_mapping.size());
    for(int i = 0; i < order_of_species_mapping.size(); i++ ) {
        ASSERT_EQ(final_order_of_species_mapping[i], order_of_species_mapping[i]);
    }
    ASSERT_EQ(finalStringStartPositions.size(), next_gene_locations.size());
    for(int i = 0; i < next_gene_locations.size(); i++ ) {
        ASSERT_EQ(finalStringStartPositions[i], next_gene_locations[i]);
    }
}

TEST_F (MotifIteratorTest, IteratoreNoCount) { // make sure print out is not binary!
    int type = 1;
    Alphabet alphabet = (Alphabet)0;
    std::ostringstream stream;
    int count = ST->printMotifs(l, alphabet, 0, *bls, stream, type == 0); // 0 == AB, 1 is AF
    std::string output = stream.str();
    int start = 0;
    int end = output.find_first_of('\n', start);
    while(end != std::string::npos) {
        std::string line = output.substr(start, (end - start));
        // verify line!
        int motifstart = line.find_first_of('\t', 0);
        int motifend = line.find_first_of('\t', motifstart + 1);
        std::string motif = line.substr(motifstart + 1, (motifend - motifstart) - 1);
        int blscount = std::stoi(line.substr(motifend + 1));
        occurence_bits occ = 0; // check occurence for each species+rcstring
        if(BDAndRC.find(motif) != std::string::npos)
            occ |= species_BD;
        if(OSAndRC.find(motif) != std::string::npos)
            occ |= species_OS;
        if(SBAndRC.find(motif) != std::string::npos)
            occ |= species_SB;
        if(ZMAndRC.find(motif) != std::string::npos)
            occ |= species_ZM;
        ASSERT_EQ(bls->getBLSVector(occ)[0], blscount);
        start = end + 1;
        end = output.find_first_of('\n', start);

    }
    // check for all other motifs that they are not preserved:
    for (int i =0; i < possible_motifs.size(); i++) {
        std::string motif = possible_motifs[i];
        occurence_bits occ = 0; // check occurence for each species+rcstring
        if(BDAndRC.find(motif) != std::string::npos)
            occ |= species_BD;
        if(OSAndRC.find(motif) != std::string::npos)
            occ |= species_OS;
        if(SBAndRC.find(motif) != std::string::npos)
            occ |= species_SB;
        if(ZMAndRC.find(motif) != std::string::npos)
            occ |= species_ZM;
        ASSERT_EQ(bls->getBLSVector(occ)[0], 0);
    }
}

TEST (Motif, MotifGroup) {
    std::vector<std::string> motifs{
        "ACGTACGT",
        "TTTACC",
        "GGTAAA",
        "GTACGTAK" // RC is MTACGTAC -> grp is AACCGMTT
    };
    std::vector<std::string> groups{
        "AACCGGTT",
        "ACCTTT",
        "AAAGGT",
        "AACGGKTT" // rc = AAMCCGTT -> group is AACCGMMT
    };
    std::vector<bool> isGroupRep{
        true,
        false,
        true,
        false,
    };
    for(int i = 0; i < motifs.size(); i++) {
        std::string grp = Motif::getGroupID(motifs[i]);
        ASSERT_EQ(Motif::isGroupRepresentative(motifs[i]), isGroupRep[i]);

    }
}

TEST_F (MotifIteratorTest, IteratoreNoCountDegenerate) { // make sure print out is not binary!
    int type = 1;
    Alphabet alphabet = (Alphabet)2;
    std::ostringstream stream;
    int count = ST->printMotifs(l, alphabet, 1, *bls, stream, type == 0); // 0 == AB, 1 is AF
    std::string output = stream.str();
    int start = 0;
    int end = output.find_first_of('\n', start);
    std::vector<std::string> foundmotifs;
    while(end != std::string::npos) {
        std::string line = output.substr(start, (end - start));
        // verify line!
        int motifstart = line.find_first_of('\t', 0);
        int motifend = line.find_first_of('\t', motifstart + 1);
        std::string motif = line.substr(motifstart + 1, (motifend - motifstart) - 1);
        ASSERT_EQ(Motif::isGroupRepresentative(motif), true);
        foundmotifs.push_back(motif);
        int blscount = std::stoi(line.substr(motifend + 1));
        std::vector<std::string> expandedmotifs = expand(motif);
        occurence_bits occ = 0; // check occurence for each species+rcstring
        for (int i = 0; i < expandedmotifs.size(); i++) {
            // std::cout << "\t" << expandedmotifs[i] << std::endl;
            if(BDAndRC.find(expandedmotifs[i]) != std::string::npos)
                occ |= species_BD;
            if(OSAndRC.find(expandedmotifs[i]) != std::string::npos)
                occ |= species_OS;
            if(SBAndRC.find(expandedmotifs[i]) != std::string::npos)
                occ |= species_SB;
            if(ZMAndRC.find(expandedmotifs[i]) != std::string::npos)
                occ |= species_ZM;
        }
        ASSERT_EQ(bls->getBLSVector(occ)[0], blscount);
        start = end + 1;
        end = output.find_first_of('\n', start);
    }
    // check for all other motifs that they are not preserved:
    for (int i =0; i < possible_motifs.size(); i++) { // for every motif i
        for(int j = 0; j < possible_motifs[i].length(); j++) { // every position j
            for(int k = 0; k < twofoldandN.length(); k++) { // replace with a twofoldandN
                std::string newmotif = possible_motifs[i];
                newmotif[j] = twofoldandN[k];
                std::vector<std::string> expandedmotifs = expand(newmotif);
                occurence_bits occ = 0;
                for (int e = 0; e < expandedmotifs.size(); e++) {
                    if(BDAndRC.find(expandedmotifs[e]) != std::string::npos)
                        occ |= species_BD;
                    if(OSAndRC.find(expandedmotifs[e]) != std::string::npos)
                        occ |= species_OS;
                    if(SBAndRC.find(expandedmotifs[e]) != std::string::npos)
                        occ |= species_SB;
                    if(ZMAndRC.find(expandedmotifs[e]) != std::string::npos)
                        occ |= species_ZM;
                }
                int blscount = bls->getBLSVector(occ)[0];
                if(blscount > 0) {
                    ASSERT_EQ(std::find(foundmotifs.begin(), foundmotifs.end(), newmotif) == foundmotifs.end() && Motif::isGroupRepresentative(newmotif), false);
                }
            }
        }
    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

#include "gtest/gtest.h"
#include <bits/stdc++.h>
#include <vector>
#include "motif.h"


class BlsTreeScoresTest: public ::testing::Test {
protected:
    int N = 4;
    std::vector<float> blsThresholds;
    occurence_bits species_BD = 1;
    occurence_bits species_OS = 2;
    occurence_bits species_SB = 4;
    occurence_bits species_ZM = 8;
    // newick without paralogue -> this should be default since I merge genes in the string
    BLSScore *bls = NULL;
    std::vector<std::string> order_of_species;
    std::vector<occurence_bits> all_possible_occurence_bits;
    std::vector<float> all_possible_scores;
    std::vector<bool> all_possible_gt_min_threshold;
    std::vector<char> all_possible_print_blsvector_binary;
    std::vector<char> all_possible_print_blsvector;
    std::vector<std::vector<bool>> all_possible_gt_threshold_i;
    std::vector<std::vector<blscounttype>> all_possible_bls_vectors;
    std::vector<std::vector<blscounttype>> cum_bls_vectors;
    std::string newick = "((BD:0.2688,OS:0.2688):0.0538,(SB:0.086,ZM:0.086):0.2366);";

    BlsTreeScoresTest( ) {
    }

    void SetUp( ) {
        // code here will execute just before the test ensues
        // initialization code here
        blsThresholds.push_back(0.15);
        blsThresholds.push_back(0.5);
        blsThresholds.push_back(0.6);
        blsThresholds.push_back(0.7);
        blsThresholds.push_back(0.9);
        blsThresholds.push_back(0.95);

        all_possible_occurence_bits.push_back(species_BD);
        all_possible_scores.push_back(0.0f);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{false, false, false, false, false, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{0, 0, 0, 0, 0, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{0, 0, 0, 0, 0, 0});
        all_possible_print_blsvector_binary.push_back(0);
        all_possible_print_blsvector.push_back('0');
        all_possible_gt_min_threshold.push_back(false);
        all_possible_occurence_bits.push_back(species_OS);
        all_possible_scores.push_back(0.0f);
        all_possible_gt_min_threshold.push_back(false);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{false, false, false, false, false, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{0, 0, 0, 0, 0, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{0, 0, 0, 0, 0, 0});
        all_possible_print_blsvector_binary.push_back(0);
        all_possible_print_blsvector.push_back('0');
        all_possible_occurence_bits.push_back(species_SB);
        all_possible_scores.push_back(0.0f);
        all_possible_gt_min_threshold.push_back(false);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{false, false, false, false, false, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{0, 0, 0, 0, 0, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{0, 0, 0, 0, 0, 0});
        all_possible_print_blsvector_binary.push_back(0);
        all_possible_print_blsvector.push_back('0');
        all_possible_occurence_bits.push_back(species_ZM);
        all_possible_scores.push_back(0.0f);
        all_possible_gt_min_threshold.push_back(false);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{false, false, false, false, false, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{0, 0, 0, 0, 0, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{0, 0, 0, 0, 0, 0});
        all_possible_print_blsvector_binary.push_back(0);
        all_possible_print_blsvector.push_back('0');
        all_possible_occurence_bits.push_back(species_BD | species_OS);
        all_possible_scores.push_back(0.5376f);
        all_possible_gt_min_threshold.push_back(true);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{true, true, false, false, false, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{1, 1, 0, 0, 0, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{1, 1, 0, 0, 0, 0});
        all_possible_print_blsvector_binary.push_back(2);
        all_possible_print_blsvector.push_back('2');
        all_possible_occurence_bits.push_back(species_BD | species_SB);
        all_possible_scores.push_back(0.6452);
        all_possible_gt_min_threshold.push_back(true);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{true, true, true, false, false, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{1, 1, 1, 0, 0, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{2, 2, 1, 0, 0, 0});
        all_possible_print_blsvector_binary.push_back(3);
        all_possible_print_blsvector.push_back('3');
        all_possible_occurence_bits.push_back(species_BD | species_ZM);
        all_possible_scores.push_back(0.6452);
        all_possible_gt_min_threshold.push_back(true);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{true, true, true, false, false, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{1, 1, 1, 0, 0, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{3, 3, 2, 0, 0, 0});
        all_possible_print_blsvector_binary.push_back(3);
        all_possible_print_blsvector.push_back('3');
        all_possible_occurence_bits.push_back(species_OS | species_SB);
        all_possible_scores.push_back(0.6452);
        all_possible_gt_min_threshold.push_back(true);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{true, true, true, false, false, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{1, 1, 1, 0, 0, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{4, 4, 3, 0, 0, 0});
        all_possible_print_blsvector_binary.push_back(3);
        all_possible_print_blsvector.push_back('3');
        all_possible_occurence_bits.push_back(species_OS | species_ZM);
        all_possible_scores.push_back(0.6452);
        all_possible_gt_min_threshold.push_back(true);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{true, true, true, false, false, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{1, 1, 1, 0, 0, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{5, 5, 4, 0, 0, 0});
        all_possible_print_blsvector_binary.push_back(3);
        all_possible_print_blsvector.push_back('3');
        all_possible_occurence_bits.push_back(species_SB | species_ZM);
        all_possible_scores.push_back(0.172f);
        all_possible_gt_min_threshold.push_back(true);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{true, false, false, false, false, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{1, 0, 0, 0, 0, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{6, 5, 4, 0, 0, 0});
        all_possible_print_blsvector_binary.push_back(1);
        all_possible_print_blsvector.push_back('1');
        all_possible_occurence_bits.push_back(species_BD | species_OS | species_SB);
        all_possible_scores.push_back(0.914f);
        all_possible_gt_min_threshold.push_back(true);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{true, true, true, true, true, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{1, 1, 1, 1, 1, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{7, 6, 5, 1, 1, 0});
        all_possible_print_blsvector_binary.push_back(5);
        all_possible_print_blsvector.push_back('5');
        all_possible_occurence_bits.push_back(species_BD | species_OS | species_ZM);
        all_possible_scores.push_back(0.914f);
        all_possible_gt_min_threshold.push_back(true);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{true, true, true, true, true, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{1, 1, 1, 1, 1, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{8, 7, 6, 2, 2, 0});
        all_possible_print_blsvector_binary.push_back(5);
        all_possible_print_blsvector.push_back('5');
        all_possible_occurence_bits.push_back(species_BD | species_SB | species_ZM);
        all_possible_scores.push_back(0.7312f);
        all_possible_gt_min_threshold.push_back(true);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{true, true, true, true, false, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{1, 1, 1, 1, 0, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{9, 8, 7, 3, 2, 0});
        all_possible_print_blsvector_binary.push_back(4);
        all_possible_print_blsvector.push_back('4');
        all_possible_occurence_bits.push_back(species_OS | species_SB | species_ZM);
        all_possible_scores.push_back(0.7312f);
        all_possible_gt_min_threshold.push_back(true);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{true, true, true, true, false, false});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{1, 1, 1, 1, 0, 0});
        cum_bls_vectors.push_back(std::vector<blscounttype>{10, 9, 8, 4, 2, 0});
        all_possible_print_blsvector_binary.push_back(4);
        all_possible_print_blsvector.push_back('4');
        all_possible_occurence_bits.push_back(species_BD | species_OS | species_SB | species_ZM);
        all_possible_scores.push_back(1.0f);
        all_possible_gt_min_threshold.push_back(true);
        all_possible_gt_threshold_i.push_back(std::vector<bool>{true, true, true, true, true, true});
        all_possible_bls_vectors.push_back(std::vector<blscounttype>{1, 1, 1, 1, 1, 1});
        cum_bls_vectors.push_back(std::vector<blscounttype>{11, 10, 9, 5, 3, 1});
        all_possible_print_blsvector_binary.push_back(6);
        all_possible_print_blsvector.push_back('6');

        bls = new BLSScore(blsThresholds, newick, N, order_of_species);
    }

    void TearDown( ) {
        // code here will be called just after the test completes
        // ok to through exceptions from here if need be
        delete bls;
    }

    ~BlsTreeScoresTest( )  {
        // cleanup any pending stuff, but no exceptions allowed
    }
    // put in any custom data members that you need
};

TEST_F (BlsTreeScoresTest, OrderOfSpecies) {
    ASSERT_EQ (order_of_species.size() , 4);
    ASSERT_EQ (order_of_species[0] , "BD");
    ASSERT_EQ (order_of_species[1] , "OS");
    ASSERT_EQ (order_of_species[2] , "SB");
    ASSERT_EQ (order_of_species[3] , "ZM");
    // check bls vector size:
    ASSERT_EQ (bls->getBLSVectorSize(), 6);
}

TEST_F (BlsTreeScoresTest, BlsScores) {
    blscounttype *cumvector = NULL;
    for (int i = 0; i < all_possible_scores.size(); i++) {
        ASSERT_FLOAT_EQ(bls->getBLSScore(all_possible_occurence_bits[i]), all_possible_scores[i]);
        ASSERT_FLOAT_EQ(bls->getBLSVector(all_possible_occurence_bits[i])[0], all_possible_print_blsvector_binary[i]);
        ASSERT_EQ(bls->greaterThanMinThreshold(all_possible_occurence_bits[i]), all_possible_gt_min_threshold[i]);
        blscounttype *vector = bls->createBlsVectorFromByte(all_possible_occurence_bits[i]);
        if(i == 0)
            cumvector = vector;
        else
            bls->addByteToBlsVector(cumvector, all_possible_occurence_bits[i]);
        for (int j = 0; j < blsThresholds.size(); j++) {
            ASSERT_EQ(bls->greaterThanThreshold(all_possible_occurence_bits[i], j), all_possible_gt_threshold_i[i][j]);
            ASSERT_EQ(vector[j], all_possible_bls_vectors[i][j]);
            ASSERT_EQ(cumvector[j], cum_bls_vectors[i][j]);
        }
        if(i>0)
            delete vector;

    }
    delete cumvector;
}

TEST_F (BlsTreeScoresTest, PrintToStream) {
    for (int i = 0; i < all_possible_scores.size(); i++) {
        std::ostringstream stream;
        std::ostringstream binstream;
        bls->writeBLSVectorInBinary(all_possible_occurence_bits[i], binstream);
        bls->writeBLSVector(all_possible_occurence_bits[i], stream);
        ASSERT_EQ( binstream.str().c_str()[0], all_possible_print_blsvector_binary[i]);
        ASSERT_EQ( stream.str().c_str()[0], all_possible_print_blsvector[i]);
    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

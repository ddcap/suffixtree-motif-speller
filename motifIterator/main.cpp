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
#include <fstream>
#include <bitset>
#include "suffixtree.h"
#include "genefamily.h"

using namespace std;

int main(int argc, char* argv[])
{
        if (argc == 8 || argc == 9) {
            int mode = 0; // motif discovery
            int type = -1; // error if not given properly!
            if (strcmp(argv[2], "AB") == 0) { type = 0; std::cerr << "Alignment Based" << std::endl; };
            if (strcmp(argv[2], "AF") == 0) { type = 1; std::cerr << "Alignment Free" << std::endl; };
            // if (strcmp(argv[2], "-") == 0) { type = 2; std::cerr << "Find Motifs" << std::endl; };
            Alphabet alphabet = (Alphabet)std::stoi(argv[3]);
            std::vector<float> blsThresholds;
            std::string s = argv[4];
            std::string delimiter = ",";

            size_t pos = 0;
            std::string token;
            std::cerr << "BLS thresholds: ";
            while ((pos = s.find(delimiter)) != std::string::npos) {
                token = s.substr(0, pos);
                blsThresholds.push_back(std::stof(token));
                std::cerr << token << ", ";
                s.erase(0, pos + delimiter.length());
            }
            blsThresholds.push_back(std::stof(s));
            std::cerr << s << std::endl;

            int maxDegeneration = std::stoi(argv[5]);
            std::pair<short, short> l(std::stoi(argv[6]), std::stoi(argv[7]));
            assert(l.first < l.second); // else empty range!
            bool countBls = (argc == 9 ? (strcmp(argv[8], "true") == 0) : false);

            if ((strcmp(argv[1], "-") == 0))
                GeneFamily::readOrthologousFamily(mode, std::cin, blsThresholds, alphabet, type, l, maxDegeneration, countBls);
            else
                GeneFamily::readOrthologousFamily(mode, argv[1], blsThresholds, alphabet, type, l, maxDegeneration, countBls);
        } else if (argc == 6 || argc == 7) {
            int mode = 1; // find motif location
            int type = -1; // error if not given properly!
            if (strcmp(argv[2], "AB") == 0) { type = 0; std::cerr << "Alignment Based" << std::endl; };
            if (strcmp(argv[2], "AF") == 0) { type = 1; std::cerr << "Alignment Free" << std::endl; };
            Alphabet alphabet = (Alphabet)3;
            std::vector<float> blsThresholds;
            std::string s = argv[3];
            std::string delimiter = ",";

            size_t pos = 0;
            std::string token;
            std::cerr << "BLS thresholds: ";
            while ((pos = s.find(delimiter)) != std::string::npos) {
                token = s.substr(0, pos);
                blsThresholds.push_back(std::stof(token));
                std::cerr << token << ", ";
                s.erase(0, pos + delimiter.length());
            }
            blsThresholds.push_back(std::stof(s));
            std::cerr << s << std::endl;

            int maxDegeneration = std::stoi(argv[4]);
            int maxLen = std::stoi(argv[5]);
            std::pair<short, short> l(-1, maxLen);
            float min_bls = (argc == 7 ? std::stof(argv[6]) : 0);

            if ((strcmp(argv[1], "-") == 0))
                GeneFamily::readOrthologousFamily(mode, std::cin, blsThresholds, alphabet, type, l, maxDegeneration, false, min_bls);
            else
                GeneFamily::readOrthologousFamily(mode, argv[1], blsThresholds, alphabet, type, l, maxDegeneration, false, min_bls);
        } else {
            std::cerr << "usage: " << std::endl;
            std::cerr << "DISCOVERY: ./motifIterator input type alphabet blsThresholdList degeneration minlen maxlen [countBls]" << std::endl;
            std::cerr << "\tinput:\tInput file or '-' for stdin." << std::endl;
            std::cerr << "\ttype:\tAB or AF for alignment based or alignment free motif discovery" << std::endl;
            std::cerr << "\talphabet (int):\t0: Exact, 1: Exact And N, 2: Exact, Twofolds And N, 3: All" << std::endl;
            std::cerr << "\tblsThresholdList:\tComma sepparated list of bls thresholds (between 0 to 1). Example '0.15,0.5,0.6,0.7,0.9,0.95'" << std::endl;
            std::cerr << "\tdegeneration:\tNumber of degenerate characters." << std::endl;
            std::cerr << "\tminlen:\tMinimum motif length, inclusive (i.e. length >= minlen)." << std::endl;
            std::cerr << "\tmaxlen:\tMaximum motif length, non inclusive (i.e. length < maxlen)." << std::endl;
            std::cerr << "\tcountBls:\tIndicates whether valid motifs per BLS threshold should be counted. true or [false]." << std::endl;
            std::cerr << "MATCH MOTIFS: ./motifIterator input type blsThresholdList degeneration maxlen [bls_threshold]" << std::endl;
            std::cerr << "\tinput:\tInput file or '-' for stdin: ortho group file followed by a list of sorted motifs to find" << std::endl;
            std::cerr << "\ttype:\tAB or AF for alignment based or alignment free motif discovery" << std::endl;
            std::cerr << "\tblsThresholdList:\tComma sepparated list of bls thresholds (between 0 to 1). Example '0.15,0.5,0.6,0.7,0.9,0.95'" << std::endl;
            std::cerr << "\tdegeneration:\tNumber of degenerate characters." << std::endl;
            std::cerr << "\tmaxlen:\tMaximum motif length, non inclusive (i.e. length < maxlen)." << std::endl;
            return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
}

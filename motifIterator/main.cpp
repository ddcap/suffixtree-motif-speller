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

        if(argc < 7) {
            std::cerr << "usage: ./motifIterator input type alphabet blsThresholdList degeneration minlen maxlen" << std::endl;
            std::cerr << "\tinput:\tInput file or '-' for stdin." << std::endl;
            std::cerr << "\ttype:\tAB or AF" << std::endl;
            std::cerr << "\talphabet (int):\t0: Exact, 1: Exact And N, 2: Exact, Twofolds And N, 3: All" << std::endl;
            std::cerr << "\tblsThrehsoldList:\tComma sepparated list of bls thresholds (between 0 to 1). Example '0.15,0.5,0.6,0.7,0.9,0.95'" << std::endl;
            std::cerr << "\tdegeneration:\tNumber of degenerate characters." << std::endl;
            std::cerr << "\tminlen:\tMinimum motif length, inclusive (i.e. length >= minlen)." << std::endl;
            std::cerr << "\tmaxlen:\tMaximum motif length, non inclusive (i.e. length < maxlen)." << std::endl;
            return EXIT_FAILURE;
        }
        bool typeIsAB = (strcmp(argv[2], "AB") == 0);
        std::cerr << (typeIsAB ? "Alignment Based" : "Alignment Free") << std::endl;
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

        if ((strcmp(argv[1], "-") == 0))
            GeneFamily::readOrthologousFamily(std::cin, blsThresholds, alphabet, typeIsAB, l, maxDegeneration);
        else
            GeneFamily::readOrthologousFamily(argv[1], blsThresholds, alphabet, typeIsAB, l, maxDegeneration);
        return EXIT_SUCCESS;
}
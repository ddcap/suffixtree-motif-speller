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
        std::pair<short, short> l(6, 13);
        int maxDegeneration = 3;
        if(argc < 3) {
            std::cerr << "usage: ./motifIterator type alphabet [file/stdin] [maxDegeneration] [l.min l.max]" << std::endl;
            std::cerr << "\ttype: AB or AF" << std::endl;
            std::cerr << "\talpbat (int):\t0: Exact, 1: Exact And N, 2: Exact, Twofolds And N, 3: All" << std::endl;
        }
        bool typeIsAB = (strcmp(argv[1], "AB") == 0);
        Alphabet alphabet = (Alphabet)std::stoi(argv[2]);
        std::cerr << (typeIsAB ? "Alignment Based" : "Alignment Free") << std::endl;
        if (argc >= 3) {
            if (argc >= 4 ) {
                maxDegeneration = std::stoi(argv[4]);
                if (argc >= 6) {
                    l.first = std::stoi(argv[5]);
                    l.second = std::stoi(argv[6]);
                }
            }

            if ((strcmp(argv[3], "stdin") == 0) || (strcmp(argv[3], "-") == 0))
                GeneFamily::readOrthologousFamily(std::cin, alphabet, typeIsAB, l, maxDegeneration);
            else
                GeneFamily::readOrthologousFamily(argv[3], alphabet, typeIsAB, l, maxDegeneration);
        } else {
            GeneFamily::readOrthologousFamily(std::cin, alphabet, typeIsAB, l, maxDegeneration);
        }
        return EXIT_SUCCESS;
}

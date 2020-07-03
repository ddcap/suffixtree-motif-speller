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
#include "genefamily.h"

using namespace std;


void readInput(const string& filename, string& T)
{
        ifstream ifs(filename.c_str());

        while (ifs) {
                string line;
                getline(ifs, line);

                if (line.empty())
                        continue;

                if (!T.empty())
                        T.push_back('#');
                T.append(line);
        }

        T.push_back('$');
}

int main(int argc, char* argv[])
{
        std::pair<short, short> l(6, 13);
        int maxDegeneration = 3;
        if (argc >= 2) {
            if (argc >= 3 ) {
                maxDegeneration = std::stoi(argv[2]);
                if (argc >= 5) {
                    l.first = std::stoi(argv[3]);
                    l.second =std::stoi(argv[4]);
                }
            }
            GeneFamily::readOrthologousFamily(argv[1], l, maxDegeneration);
        } else {
            GeneFamily::readOrthologousFamily(std::cin, l, maxDegeneration);
        }
        return EXIT_SUCCESS;
}

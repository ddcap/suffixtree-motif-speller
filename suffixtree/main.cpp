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
        if (argc == 2) {
            GeneFamily::readOrthologousFamily(argv[1]);
        } else {
            GeneFamily::readOrthologousFamily(std::cin);
        }
        return EXIT_SUCCESS;
}

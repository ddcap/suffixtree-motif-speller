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
#include <unordered_set>
#include "newick.h"
#include "genes.h"
#include "orthology.h"

using namespace std;

int main(int argc, char* argv[])
{
    if(argc != 4) {
        std::cerr << "usage: ./prepInput [folder of fasta files] [orthology file] [newick lengths file]" << std::endl;
        std::cerr << argv[0] << std::endl;
    }

    std::cerr << "reading fasta files from '" << argv[1] << "'" << std::endl;
    Genes genes(argv[1]);
    std::cerr << "reading tree branch lengths from '" << argv[3] << "'" << std::endl;
    Newick newick(argv[3]);
    std::cerr << "reading orthology from '" << argv[2] << "'" << std::endl;
    Orthology orthology(&genes, &newick, argv[2]);

    return EXIT_SUCCESS;
}

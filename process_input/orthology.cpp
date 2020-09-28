
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <unordered_map>
#include "orthology.h"
#include "genes.h"
#include "newick.h"

void Orthology::formatGeneList(std::ostream& o, Genes *genemap, Newick *newick, std::string cluster, std::string genelist) {
    std::cerr << "cluster " << cluster << std::endl;
    o << cluster << '\n';
    std::unordered_map<std::string, std::string> species_to_genes_map;
    std::unordered_map<std::string, std::string> species_to_geneids_map;
    std::istringstream iss(genelist);
    std::string item;
    while (std::getline(iss, item, delim)) {
        // std::cerr << "item " << item;
        Gene *gene = genemap->getGene(item);
        // std::cerr << " gene: " << gene->species << std::endl;
        if(gene == NULL) { std::cerr << "gene " + item + " not found"  << std::endl; }
        else {
            auto findpos = species_to_genes_map.find(gene->species);
            if(findpos == species_to_genes_map.end()) {
                species_to_genes_map[gene->species] = gene->sequence;
                species_to_geneids_map[gene->species] = item;
            } else {
                species_to_genes_map[gene->species] += delim + gene->sequence;
                species_to_geneids_map[gene->species] += delim + item;
            }
        }
    }
    std::vector<std::string> keys;
    keys.reserve(species_to_genes_map.size());
    for (auto x : species_to_genes_map) {
        // std::cerr << x.first << std::endl;
        keys.push_back(x.first);
    }
    newick->print_renormalised_tree(o, keys);
    o << keys.size() << '\n';
    for (auto x : species_to_genes_map) {
        o << species_to_geneids_map[x.first] << '\t' << x.first << '\n';
        o << x.second << '\n';
    }

}
void Orthology::readOrthology(Genes *genemap, Newick *newick, std::string orthologyFile, std::string outputfolder) {

    std::ifstream f(orthologyFile);
    std::string line, cluster, genes;
    int species_count; // , gene_count;
    size_t count = 0, splitpos1, splitpos2;
    if (f.is_open()) {
        while ( getline (f, line) )  {
            if(!line.empty() && line[0] != '#') {
                splitpos1 = line.find('\t', 0);
                if(splitpos1 == std::string::npos) {
                    return;
                }
                cluster = line.substr(0, splitpos1);
                splitpos2 = line.find('\t', splitpos1+1);
                // gene_count = std::stoi(line.substr(splitpos1+1, splitpos2));
                splitpos1 = splitpos2;
                splitpos2 = line.find('\t', splitpos1+1);
                species_count = std::stoi(line.substr(splitpos1+1, splitpos2));
                if(species_count > 1) {
                    count++;
                    genes = line.substr(splitpos2+1);
                    // create ofstream o
                    std::ofstream o;
                    o.open(outputfolder + "/" + cluster);
                    formatGeneList(o, genemap, newick, cluster, genes);
                    o.close();
                }
            }
        }
    }
    f.close();
    std::cerr << "found " << count << " clusters with more than 1 species" << std::endl;
}

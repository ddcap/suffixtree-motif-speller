
#include <string>
#include <vector>
#include <fstream>
#include <vector>
#include "newick.h"

// TODO: read strings with gene id that maps to species + sequence, this way we can make a set! (unsorted_set?) of species based on the list of genes.

void Newick::print_renormalised_tree(std::ostream& o, std::vector<std::string> species) {
    // std::cerr << "printing renormalised tree with these species: " << std::endl;
    root->resetUsed();
    for (auto x : species) {
        std::cerr << "- " << x << std::endl;
        root->useSpecies(x);
        // std::cerr << "tree: \n" << *root << std::endl;
    }
    root->getMinimimTree();
    float subtree_sum = root->get_length_sum();
    // std::cerr << "sum : " << subtree_sum << std::endl;
    // std::cerr << "final tree: ";
    // root->print_renormalised(std::cerr, subtree_sum);
    // std::cerr << std::endl;
    // o << "final newick: ";
    root->print_newick_renormalised(o, subtree_sum);
    o << '\n';
    // loop over tree and connect species with minimum tree and get sum.
}
void Newick::remove_newick_nodes() {
    // Depth-first traversal of the tree
    std::stack<newick_node *> stack;
    stack.push(root);

    while (!stack.empty()) {
            newick_node* node = stack.top();
            stack.pop();
            if(node->getChild() != NULL) stack.push(node->getChild());
            if(node->getNext() != NULL) stack.push(node->getNext());

            delete node;
    }
}
void Newick::prep_newick(std::string lengthsfile) {
    std::ifstream f(lengthsfile);
    std::string line;
    if (f.is_open()) {
        while ( getline (f, line) )  {
            if(!line.empty() && line[0] != '#' ) {
                int splitpos = line.find_first_of(delim);
                std::string name = line.substr(0, splitpos);
                int len = std::stoi(line.substr(splitpos));
                sum += len;
                lengths[name] = len;
            }
        }
        f.close();
    }
    // create the tree structure, this is fixed! but lengths can differ!
    root = new newick_node("poaceae");
    newick_node* BOP_clade = root->addChild("BOP_clade", lengths["BOP_clade"]/sum);
    newick_node* oryza = BOP_clade->addChild("oryza", lengths["oryza"]/sum);
    newick_node* oryzasativa = oryza->addChild("oryzasativa", lengths["oryzasativa"]/sum);
    newick_node* osa = oryzasativa->addChild("osa", lengths["osa"]/sum);
    osa->addNext("osaindica", lengths["osaindica"]/sum, oryzasativa);
    oryzasativa->addNext("obr", lengths["obr"]/sum, oryza);
    newick_node* pooideae_and_ped = oryza->addNext("pooideae_and_ped", lengths["pooideae_and_ped"]/sum, BOP_clade);
    newick_node* pooideae = pooideae_and_ped->addChild("pooideae", lengths["pooideae"]/sum);
    pooideae->addNext("ped", lengths["ped"]/sum, pooideae_and_ped);
    newick_node* triticeae = pooideae->addChild("triticeae", lengths["triticeae"]/sum);
    triticeae->addNext("bdi", lengths["bdi"]/sum, triticeae);
    newick_node* tricium = triticeae->addChild("tricium", lengths["tricium"]/sum);
    tricium->addNext("hvu", lengths["hvu"]/sum, pooideae);
    newick_node* ttu = tricium->addChild("ttu", lengths["ttu"]/sum);
    ttu->addNext("tae", lengths["tae"]/sum, tricium);

    newick_node* PACMAD_clade = BOP_clade->addNext("PACMAD_clade", lengths["PACMAD_clade"]/sum, root);
    newick_node* chloridoideae = PACMAD_clade->addChild("chloridoideae", lengths["chloridoideae"]/sum);
    newick_node* oth = chloridoideae->addChild("oth", lengths["oth"]/sum);
    oth->addNext("zjn", lengths["zjn"]/sum, chloridoideae);
    newick_node* panicoideae = chloridoideae->addNext("panicoideae", lengths["panicoideae"]/sum, PACMAD_clade);
    newick_node* paniceae = panicoideae->addChild("paniceae", lengths["paniceae"]/sum);
    newick_node* sit = paniceae->addChild("sit", lengths["sit"]/sum);
    sit->addNext("cam", lengths["cam"]/sum, paniceae);
    newick_node* andropogoneae = paniceae->addNext("andropogoneae", lengths["andropogoneae"]/sum, panicoideae);
    newick_node* ssp = andropogoneae->addChild("ssp", lengths["ssp"]/sum);
    newick_node* zeamays = ssp->addNext("zeamays", lengths["zeamays"]/sum, andropogoneae);
    zeamays->addNext("sbi", lengths["sbi"]/sum, andropogoneae);
    newick_node* zma = zeamays->addChild("zma", lengths["zma"]/sum);
    zma->addNext("zma-ph207", lengths["zma-ph207"]/sum, zeamays);

    // std::cerr << "tree: \n" << *root << std::endl;
    // std::cerr << "newick: ";
    // root->print_newick(std::cerr);
    // std::cerr << std::endl;
    // std::cerr << "sum : " << root->get_length_sum() << std::endl;
}

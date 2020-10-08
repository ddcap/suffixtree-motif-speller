
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
        // std::cerr << "- " << x << std::endl;
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

newick_node* Newick::get_first_child(std::string &newick_string, int level) {
    size_t name_end = newick_string.find(':');
    std::string name = newick_string.substr(0, name_end);
    if(name.compare(";") == 0) name = "";
    newick_string.erase(0, name_end + 1);
    float length = 0.0f;
    // next is either a ',' or a ')'
    size_t length_end = 0;
    size_t firstcomma = newick_string.find(',');
    size_t firstbracket = newick_string.find(')');
    if(firstcomma == std::string::npos) length_end = firstbracket;
    else if(firstcomma == std::string::npos) length_end = firstcomma;
    else length_end = std::min(firstbracket, firstcomma);
    if(length_end != std::string::npos) length = std::stof(newick_string.substr(0, length_end)); // can still be npos at the end, the length should be 0!
    newick_string.erase(0, length_end);
    return new newick_node(name, length, level);
}
newick_node* Newick::build_newick_recursive(std::string &newick_string, int level) {
    // std::cerr << "l[" << level << "] new branch: " << newick_string << std::endl;
    newick_node * firstchild = NULL;
    if(newick_string[0] == '(') newick_string.erase(0, 1);
    // if(newick_string[0] == '(') {
    //     newick_string.erase(0, 1);
    //     firstchild = build_newick_recursive(newick_string, level + 1);
    // } else {
        if(newick_string[0] == '(')
            firstchild = build_newick_recursive(newick_string, level + 1);
        else
            firstchild = get_first_child(newick_string, level + 1);
        // std::cerr << "l[" << level << "] 1member of branch:\n" << *firstchild << std::endl;
        // std::cerr << "l[" << level << "] string: " << newick_string << std::endl;
        newick_node *child = firstchild;
        newick_node *next = NULL;
        while(newick_string[0] == ',') {
            newick_string.erase(0, 1);
            if(newick_string[0] == '(')
                next = build_newick_recursive(newick_string, level + 1);
            else
                next = get_first_child(newick_string, level + 1);
            // std::cerr << "l[" << level << "] 2member of branch:\n" << *next << std::endl;
            // std::cerr << "l[" << level << "] string: " << newick_string << std::endl;
            child = child->setNext(next);
        }
    // }
    // std::cerr << "l[" << level << "] children processed! " << newick_string << std::endl;
    // next should be a ')'
    if(newick_string[0] != ')') {
        std::cerr << "after child should be ) [" + newick_string + "]" << std::endl;
        throw "after child should be ) [" + newick_string + "]";
    } else {
        newick_string.erase(0, 1);
        newick_node *parent = get_first_child(newick_string, level);
        // std::cerr << "l[" << level << "] current branch node:\n" << *parent << std::endl;
        // std::cerr << "updating child parent " << *firstchild << std::endl;
        parent->setChild(firstchild);
        firstchild->setParent(parent);
        while(firstchild->getNext() != NULL) {
            firstchild = firstchild->getNext();
            // std::cerr << "updating child parent " << *firstchild << std::endl;
            firstchild->setParent(parent);
        }
        // update children with parent node and level!
        // std::cerr << "l[" << level << "] current branch node:\n" << *parent << std::endl;
        return parent;
    }
}
void Newick::prep_newick(std::string newickfile) {
    std::ifstream f(newickfile);
    std::string line;
    if (f.is_open()) {
        while ( root == NULL && getline (f, line) )  {
            if(!line.empty() && line[0] != '#' ) {
                try{
                    root = build_newick_recursive(line, 0);
                } catch (const char* msg) {
                    std::cerr << "[ERROR] " << msg << std::endl;
                }
            }
        }
        f.close();
    }
    // std::cerr << "tree: \n" << *root << std::endl;
    // std::vector<std::string> keys = {
    //     "zma", "zma-ph207", "bdi", "cam",
    //     "hvu", "obr", "osa", "osaindica",
    //     "oth", "ped", "sbi", "sit",
    //     "ssp", "tae", "ttu", "zjn"};
    // print_renormalised_tree(std::cerr, keys);
    // std::cerr << "newick: ";
    // root->print_newick(std::cerr);
    // std::cerr << std::endl;
    // std::cerr << "sum : " << root->get_length_sum() << std::endl;
}

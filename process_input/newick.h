#ifndef NEWICK_H
#define NEWICK_H


#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <list>
#include <stack>

class newick_node {
private:
    std::string name;
    float length;
    int level;
    bool used;
    newick_node* parent;
    newick_node* next;
    newick_node* child;
    std::ostream& write(std::ostream& o) const{
        if(used) {
            for(int i = 0; i < level; i++ ) { o << " -"; }
            if(level > 0) o << name << " [" << length << "]"; // print this actual node
            else o << name;
        }
        if(child != NULL) {
            if(child->used) o << std::endl;
            o << *child;
        }
        if(next != NULL) {
            if(next->used) o << std::endl;
            o << *next;
        }
        return o;
    }
    std::ostream&  write_renormalised(std::ostream& o, float sum) const{
        if(used) {
            for(int i = 0; i < level; i++ ) { o << " -"; }
            if(level > 0) o << name << " [" << length /sum<< "]"; // print this actual node
            else o << name;
        }
        if(child != NULL) {
            if(child->used) o << std::endl;
            child->write_renormalised(o, sum);
        }
        if(next != NULL) {
            if(next->used) o << std::endl;
            next->write_renormalised(o, sum);
        }
        return o;
    }
    std::ostream&  write_newick(std::ostream& o) const{
        // if(!used) return o;
        if(child != NULL) {
            if(used) o << "(";
            child->write_newick(o);
            if(used) o << ")";
        } else {
            if(used) o << name; // print this actual node
        }
        if(used && level != 0) o << ":" << length; // print length
        if(next != NULL) {
            if(used && next->used) o << ",";
            next->write_newick(o);
        }
        return o;
    }
    std::ostream&  write_newick_renormalised(std::ostream& o, float sum) const{
        // if(!used) return o;
        if(child != NULL) {
            if(used) o << "(";
            child->write_newick_renormalised(o, sum);
            if(used) o << ")";
        } else {
            if(used) o << name; // print this actual node
        }
        if(used && level != 0) o << ":" << length / sum; // print length
        if(next != NULL) {
            if(used && next->used) o << ",";
            next->write_newick_renormalised(o, sum);
        }
        return o;
    }
    newick_node(std::string name_, float length_, int level_, newick_node *parent_) : name(name_), length(length_), level(level_), used(true), parent(parent_), next(NULL), child(NULL) {
    }
    newick_node *recFindSpecies(std::string x) {
        if(x.compare(name) == 0)
            return this;
        else {
            if(child != NULL) {
                newick_node* node = child->recFindSpecies(x);
                if (node != NULL) return node;
            }
            if(next != NULL) {
                newick_node* node = next->recFindSpecies(x);
                if (node != NULL) return node;
            }
            return NULL; // in case its not found!
        }
    }
public:
    newick_node(std::string name_) : name(name_), length(0), level(0), used(false), parent(NULL), next(NULL), child(NULL) {
    }
    newick_node *addNext(std::string name_, float length_, newick_node *parent_) { next = new newick_node(name_, length_, level, parent_); return next;}
    newick_node *addChild(std::string name_, float length_) { child = new newick_node(name_, length_, level + 1, this); return child; }
    newick_node *getChild() const { return child; }
    newick_node *getNext() const { return next; }
    bool useSpecies(std::string x) {
        // find the species, depth first.
        newick_node *found_node = recFindSpecies(x);
        if(found_node == NULL) {
            std::cerr << "species " << x  << " not found, skipping " << std::endl;
            return false;
        }
        found_node->used = true;
        newick_node *found_node_parent = found_node->getParent();
        // go up to root again and check if siblings are used, if siblings are used, parent must be set to used as well!
        while(found_node_parent != NULL) {
            found_node_parent->used = true;
            found_node_parent = found_node_parent->getParent();
        }
        return true;
    }
    void getMinimimTree() {
        int siblingsUsed = 0;
        if(used) siblingsUsed += 1;
        newick_node* nextSibling = next;
        while(nextSibling != NULL) {
            if(nextSibling->used) siblingsUsed += 1;
            nextSibling = nextSibling->next;
        }
        // std::cerr << name << " " << siblingsUsed << " used siblings " << std::endl;
        if (siblingsUsed < 2) {
            used = false;
            // std::cerr << name << " false " << std::endl;
            newick_node* nextSibling = next;
            if(child != NULL) {
                child->getMinimimTree();
            }
            while(nextSibling != NULL) {
                nextSibling->used = false;
                if(nextSibling->child != NULL) {
                    nextSibling->child->getMinimimTree();
                }
                nextSibling = nextSibling->next;
            }
        // } else {
            // std::cerr << name << " true " << std::endl;
        }
    }
    void resetUsed() {
        used = false;
        if(child != NULL) child->resetUsed();
        if(next != NULL) next->resetUsed();
    }
    newick_node *getParent() const { return parent; }
    float get_length_sum() {
        float sum = 0.0f;
        if(used) sum += length;
        if(child!=NULL) sum += child->get_length_sum();
        if(next!=NULL) sum += next->get_length_sum();
        return sum;
    }

    void print_newick(std::ostream& o) {
        o << "(";
        this->write_newick(o);
        o << ");";
    }
    void print_newick_renormalised(std::ostream& o, float sum) {
        o << "(";
        this->write_newick_renormalised(o, sum);
        o << ");";
    }
    friend std::ostream& operator<< (std::ostream& o, const newick_node& b) {
        return b.write(o);
    }
    void print_renormalised(std::ostream& o, float sum) {
        this->write_renormalised(o, sum);
    }
};

class Newick {
private:
  float sum; // this should be adapted to the sum required to connect all branches that are present...
  std::unordered_map<std::string, int> lengths;
  char delim = '\t';
  newick_node* root;

  void prep_newick(std::string lengthsfile);
  void remove_newick_nodes();
public:
    Newick(std::string lengthsfile) : sum(0){
        prep_newick(lengthsfile);
    }
    ~Newick() {
        remove_newick_nodes();
    }
    void print_renormalised_tree(std::ostream& o, std::vector<std::string> species);
};


#endif

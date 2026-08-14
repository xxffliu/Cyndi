#include "../include/Fragment.h"

// constructors
FRAGMENT::FRAGMENT():
                     id_(0),
                     name_("UNK"),
                     num_mem_atoms_(0),
                     num_mem_bonds_(0),
                     num_root_node_atoms_(0),
                     root_node_atom_(),
                     memberatom_(),
                     memberbond_(){}
                     
FRAGMENT::FRAGMENT(const FRAGMENT& frag):
                         id_(frag.id_),
                         name_(frag.name_),
                         num_mem_atoms_(frag.num_mem_atoms_),
                         num_mem_bonds_(frag.num_mem_bonds_),
                         num_root_node_atoms_(frag.num_root_node_atoms_),
                         root_node_atom_(frag.root_node_atom_),
                         memberatom_(frag.memberatom_),
                         memberbond_(frag.memberbond_){}
                         
FRAGMENT& FRAGMENT::operator=(const FRAGMENT& frag){
          id_ = frag.id_;
          name_ = frag.name_;
          num_mem_atoms_ = frag.num_mem_atoms_;
          num_mem_bonds_ = frag.num_mem_bonds_;
          num_root_node_atoms_ = frag.num_root_node_atoms_;
          root_node_atom_ = frag.root_node_atom_;
          memberatom_ = frag.memberatom_;
          memberbond_ = frag.memberbond_;
		  return *this;
          }

void FRAGMENT::clear(){
     id_ = 0;
     name_ = "UNK";
     num_mem_atoms_ = 0;
     num_mem_bonds_ = 0;
     num_root_node_atoms_ = 0;
     root_node_atom_.clear();
     memberatom_.clear();
     memberbond_.clear();
     }
     
 FRAGMENT::~FRAGMENT(){
                       clear();
 }
 
 // get and set methods
int FRAGMENT::get_id() const{
    return id_;
}

string FRAGMENT::get_name() const{
    return name_;
}

int FRAGMENT::get_num_mem_bonds() const{
    return num_mem_bonds_;
}

int FRAGMENT::get_num_mem_atoms() const{
    return num_mem_atoms_;
}

int FRAGMENT::get_num_root_node_atoms() const{
    return num_root_node_atoms_;
}

BONDVec& FRAGMENT::get_member_bonds(){
        return memberbond_;
}
             
ATOMVec& FRAGMENT::get_member_atoms(){
        return memberatom_;
}

vector<int>& FRAGMENT::get_root_node_atoms(){
    return root_node_atom_;
}                     

void FRAGMENT::set_id(const int i){
     id_ = i;
}

void FRAGMENT::set_name(const string n){
     name_ = n;
}

void FRAGMENT::set_num_root_node_atoms(const int n){
     num_root_node_atoms_ = n;
}

void FRAGMENT::set_num_mem_bonds(const int n){
     num_mem_bonds_ = n;
}

void FRAGMENT::set_num_mem_atoms(const int n){
     num_mem_atoms_ = n;
}

void FRAGMENT::set_member_bonds(const BONDVec& bvec){
     memberbond_ = bvec;
}

void FRAGMENT::set_member_atoms(const ATOMVec& avec){
     memberatom_ = avec;
}

void FRAGMENT::set_root_node_atoms(const vector<int>& atom){
     root_node_atom_ = atom;
}

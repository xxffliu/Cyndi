#ifndef Fragment_h
#define Fragment_h

#ifndef ATOM_H
#include "Atom.h"
#endif

#ifndef BOND_H
#include "Bond.h"
#endif


using namespace std;

class FRAGMENT{
      public:
             // default constructor
             FRAGMENT();
             // copy constructor
             FRAGMENT(const FRAGMENT& frag);
             // assign constructor
             FRAGMENT& operator=(const FRAGMENT& frag);
             // destructor
             virtual ~FRAGMENT();
             // a clear function
             void clear();
             
             // get and set memeber functions
             int get_id() const;
             string get_name() const;
             int get_num_mem_bonds() const;
             int get_num_mem_atoms() const;
             int get_num_root_node_atoms() const;
             BONDVec& get_member_bonds();
             ATOMVec& get_member_atoms();
             vector<int>& get_root_node_atoms();
             
             void set_id(const int i);
             void set_name(const string n);
             void set_num_root_node_atoms(const int n);
             void set_num_mem_bonds(const int n);
             void set_num_mem_atoms(const int n);
             void set_member_bonds(const BONDVec& bvec);
             void set_member_atoms(const ATOMVec& avec);
             void set_root_node_atoms(const vector<int>& atom);
             
             // 
             
      private:
              // fragment id
              int id_;
              // fagment name;
              string name_;
              // number of member atoms
              int num_mem_atoms_;
              // a vector of member atoms
              ATOMVec memberatom_;
              // number of member bonds
              int num_mem_bonds_;
              // a vector of member bonds
              BONDVec memberbond_;
              // number of root node atoms
              int num_root_node_atoms_;
              // root node atoms connecting other fragment with single rotatable bond
              vector<int> root_node_atom_;
              //string type;
              //int dict_type;
              //string chain;
              //string sub_type;
              //int inter_bond;
              //string status_bit;
};
#endif


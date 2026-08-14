#ifndef BOND_H
#define BOND_H
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <utility>

#ifndef ATOM_H
#include "Atom.h"
#endif


//using namespace std;
class ATOM;

class BOND{
      public:
             // a boolian marker for rotational bonds
             bool is_rotor;
             // constuctor list, copy constructor and assignment operator overloading
             //friend class ATOM;
			 /*// a mark flag to indicate this bond is an aromatic ring bond
			  bool  AromaticRingBond_;*/
             BOND();
             BOND(const BOND& bond0);
             BOND& operator=(const BOND& bond0);
             virtual ~BOND();
             void clear();
             // 'get' and 'set' methods
             inline ATOM* get_first_atom()
			 {
				 return first_atom;
			 }
             inline ATOM* get_second_atom()
			 {
				 return second_atom;
			 }
             ATOM* get_partner(const ATOM*);
             ATOM* get_partner(const int);
             inline int get_id() const
			 {
				 return id;
			 }
             inline string get_type() const
			 {
				 return type;
			 }
             inline float get_bond_order() const
			 {
				 return bond_order_;
			 }
             
			 bool has_atom(int id);
             //string get_status_bit();
             inline void set_first_atom(ATOM* atom0)
			 {
				 first_atom = atom0;
			 }
             inline void set_second_atom(ATOM* atom0)
			 {
				 second_atom = atom0;
			 }
             inline void set_id(int id0)
			 {
				 id = id0;
			 }
             inline void set_type(const string& type0)
			 {
				 type = type0;
			 }
             inline void set_bond_order(float order)
			 {
				 bond_order_ = order;
			 }
             //void set_status_bit(const string& status_bit0);
             // io stream
             friend istream& operator >> (istream& is, BOND& bond0);
             friend ostream& operator >> (ostream& os, const BOND& bond0);
             // other methods
             void swap(BOND& bond1,BOND& bond2);
			 inline bool isAromaticRingBond(){return AromaticRingBond_;}
             
             
             // now define some boolian function to check the bond type
             inline bool is_single()
			 {
				 return ((type == "1") && !AromaticRingBond_) ;
			 }
             inline bool is_double()
			 {
				 return (type == "2");
			 }
             inline bool is_triple()
			 {
				 return (type == "3");
			 }
             inline bool is_aromatic()
			 {
				 return (type == "ar");
			 }
             bool is_amide();
			 // check this bond is in guanidino
			 bool is_guanidino();
             //check the bond is in ring
             bool is_in_ring();
			 inline void setAromaticRingBond()
			 {
				 AromaticRingBond_ = true;
				 return ;
			 }
               
      private:
              ATOM* first_atom;
              ATOM* second_atom;
              int id;
              string type;
              float bond_order_;
			  //a boolian maker for bond which is in aromatic ring
		      bool  AromaticRingBond_;
              //string status_bit;
             
};

#endif
// define
//typedef vector<BOND*> BONDVec; 

            

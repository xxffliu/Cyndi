#ifndef ATOM_H
#define ATOM_H
#include <iostream>
#include <string>
#include <vector>
#include <utility>

#ifndef VECTOR3_H
#include "Vector3.h"
#endif

#ifndef BOND_H
#include "Bond.h"
#endif
//using namespace std;

class BOND;
class ATOM;
typedef vector<BOND*> BONDVec;
typedef vector<ATOM*> ATOMVec;
class ATOM{
	 
      public:
             friend class BOND;
			 friend class MOL2IO;
             // a boolian marker for ring atom
             bool is_ring;
             // a ring identifier
             vector<int> ring_id;
             // a boolian marker for the node atom connecting rigid fragments
             bool is_root_node;

             // a alias for any atom type
             enum{ANY_TYPE = 26};
             // a alias for maxium number of bonds
             enum{
                  MAX_NUM_BOND = 12
                  };
             inline bool operator==(const ATOM& atom)
             {
                    return name == atom.name && symbol_type_ == atom.symbol_type_
                            && _p == atom._p;
             }
             inline bool operator!=(const ATOM& atom)
             {
                    return !(*this == atom);
             }

             // constuctor list, copy constructor and assignment operator overloading
             ATOM();
             ATOM(const ATOM& atom0);
             ATOM& operator=(const ATOM& atom0);
             virtual ~ATOM();
             void clear();
             // 'get' and 'set' methods
             inline int get_id() const
			 {
				  return id;
			 }
             inline string get_name() const
			 {
				 return name;
			 }
             inline string get_symbol_type() const
			 {
				 return symbol_type_;
			 }
             inline int get_type() const
			 {
				 return type;
			 }
			 inline int get_mmff94_type() const
			 {
				 return mmff94_type;
			   
			 }
			 inline string get_mmff94_symbol_type() const
			 {
				 return mmff94_symbol_type_;
			 }
             inline int get_subst_id() const
			 {
				 return subst_id;
			 }
             inline string get_subst_name() const
			 {
				 return subst_name;
			 }
             inline float get_charge() const
			 {
				 return charge;
			 }
			 inline float get_mmff94_partial_charge()const
			 {
				 return mmff94_partial_charge;
			 }
             //int get_valence() const;
             //return the element
             inline string get_element() const
			 {
				 return element_;
			 }
             inline vector3& get_orig_position()
			 {
				 return orig_position_;
			 }
             inline vector3& get_position()
			 {
				 return _p;
			 }
             inline vector3& get_force()
			 {
				 return _f;
			 }
             // return the vdw radius of this atom
             inline float get_radius() const
			 {
				 return radius_;
			 }
			 // return the atomic weight for each element
			 inline float get_AW() const
			 {
				 return weight_;
			 }
             // get a bond vector pointer pointing to the neighbot bond list of this atom
             inline BONDVec& get_bond_list()
			 {
				 return bond_;
			 }
             // get the number of the bonds bounded to the atom
             inline int get_num_neighbor_bond() const
			 {
				 return num_of_bond_;
			 }
             // get a atom vector pointing to the neighbor atom list of this atom;
			 inline ATOMVec& get_atom_list()
			 {
				 return atom_;
			 }
             // get num of neighbor atoms
             inline int get_num_bonded_atom() const
			 {
				 return num_of_bonded_atom_;
			 }
             // get num of neighbor heavy atoms
             inline int get_num_bonded_heavy_atom() const
			 {
				 return num_of_bonded_heavy_atom_;
			 }
			 // get number of neighbor rotor bonds
			 inline int get_num_rotor_bond()
			 {
				 return num_of_rotor_bond_;
			 }
			 inline void set_num_rotor_bond(int n)
			 {
				 num_of_rotor_bond_ = n;
				 return;
			 }
			 // get the number of PI electrons
			 int get_num_PI_electrons();
             //string get_status_bit();
             inline void set_id(int id0)
			 {
				 id = id0;
			 }
             inline void set_name(const string& name0)
			 {
				 name = name0;
			 }
             //void set_coordinate(FloateVec);
             inline void set_symbol_type(const string& type0)
			 {
				  symbol_type_ = type0;
			 }
             inline void set_type(const int type0)
			 {
				 type = type0;
			 }
			 
			 inline void set_mmff94_type(const int mmff94_type0 )
			 {
				 mmff94_type=mmff94_type0;
			 
			 }
			 inline void set_mmff94_symbol_type(const string type)
			 {
				 mmff94_symbol_type_ = type;
				 return;
			 }
             inline void set_subst_id(int id0)
			 {
				 id = id0;
			 }
             inline void set_subst_name(const string& subst_name0)
			 {
				 subst_name = subst_name0;
			 }
             inline void set_charge(float charge0)
			 {
				 charge = charge0;
			 }
			 inline void set_mmff94_partial_charge(float charge0)
			 {
				 mmff94_partial_charge=charge0;
			 }
             //void set_valence(int);
             //set the element
             inline void set_element(const string& element)
			 {
				 element_ = element;
			 }
             
             inline void set_position(const vector3& pos)
			 {
				 _p = pos;
			 }
             inline void set_position(const double x,const double y,const double z)
			 {
                  _p.Set(x,y,z);
             }
             inline void set_position(const double* c)
			 {
				 _p.Set(c);
			 }
             inline void set_orig_position(vector3& pos)
			 {
				 orig_position_ =  pos;
			 }
             inline void set_force(const vector3& force)
			 {
				 _f = force;
			 }
             inline void set_force(const double x, const double y, const double z)
			 {
                 _f.Set(x,y,z);
             }
             inline void set_force(const double* f)
			 {
				 _f.Set(f);
			 }
			 // set the aromatic ring flag
			 inline void setAromaticRingAtom()
			 {
				 AromaticRingAtom_ =true;
				 return ;
			 }
             // add the bond pointer to current atom bond list
             void add_neighbor_bond_list(BOND*);
             // set the number of the bonds bounded to the atom
             void set_num_neighbor_bond(int i);
             void clear_neighbor_bond_list();
             // check the atom is bonded to the given one
             bool is_bonded_to(ATOM* atom);
             // check if the atom is geminal to the given one
             bool is_geminal_to(ATOM* atom);
             // check if the atom is vicinal to the given one
             bool is_vicinal_to(ATOM* atom);
             
             // add the neighbor bonded atoms to current atom neoghbor bonded atom list
             void add_neighbor_atom_list(ATOM*);
             // set the number of bonded atoms
             void set_num_neighbor_atom(int i);
             void set_num_neighbor_heavy_atom(int i);
             void clear_neighbor_atom_list();
             
             //void set_status_bit(const string&);
             // io stream for Mol2 file
             friend istream& operator >> (istream& is, ATOM& atom0);
             friend ostream& operator >> (ostream& os, const ATOM& atom0);
             // other methods
             void swap(ATOM& atom1, ATOM& atom2);
			                    
             // boolian methods to check the atom properities
             // a boolian function checking if this is a hydrogen atom
             inline bool is_hydrogen(){return (element_ == "H");}
             // check if is polar hydrogen
             bool is_polar_hydrogen();
             // check if is nonpolar hydrogen
             bool is_nonpolar_hydrogen();
             // check if is bonded to hydrogen atom
             bool bond_to_hydrogen();
             // check the atom types
             inline bool is_carbon(){return (element_ == "C");}
             inline bool is_nitrogen(){return (element_ == "N");}
             inline bool is_oxygen(){return (element_ == "O");}
             inline bool is_sulfur(){return (element_ == "S");}
             inline bool is_phosphorus(){return (element_ == "P");}
			 inline bool is_fluorine() {return (element_ == "F");}
			 inline bool is_clorine() {return (element_ == "Cl" || element_ == "CL");}
			 inline bool is_bromine() {return (element_ == "Br" || element_ == "BR");}
			 inline bool is_iodine() {return (element_ == "I");}
             inline bool is_aromatic(){return (symbol_type_ == "N.ar" || symbol_type_ == "C.ar");}
			 inline bool isAromaticRingAtom(){return AromaticRingAtom_;}
             inline bool is_hetero(){return (!is_carbon() && !is_hydrogen());}
             inline bool is_halogen(){return (element_ == "F" || element_ == "Cl" ||
                                       element_ == "Br" || element_ == "I");}
			 inline bool is_silicon() {return (element_ == "SI" || element_ == "Si"); }
             // check this atom is a carboxyl C (eg. CO2 or COOH)
             bool is_carboxyl();
             // check this atom is oxygen in caboxyl
             bool is_carboxyl_oxygen();
             // check this atom is a cobonyl C (eg. -CO-)
             bool is_carbonyl();
             // check if this is hydroxyl O (-OH, not COOH, SO2H, SO3H, PO3H)
             bool is_hydroxyl();
			 // check if this is sulfydryl S (-SH)
			 inline bool is_sulfydryl()
			 {
				 return (symbol_type_ == "S.3" && bond_to_hydrogen());
			 }
             // check this atom is a phosphate P (eg. PO3 or PO2)
             bool is_phosphate();
             //check this atom is oxygen in phosphate
             bool is_phosphate_oxygen();
             // check this atom is a sulfate S (eg. SO2, SO or SO3)
             bool is_sulphate();
             // check this atom is oxygen in sulfate
             bool is_sulphate_oxygen();
             // check this atom is a nitro N (eg. NO2)
             bool is_nitro();
             //check this atom is a nitro O 
             bool is_nitro_oxygen();
             // check this atom is a amide N (eg. CO-NR2)
             bool is_amide();
             // check this atom is a amine N (sp3 valence)
             bool is_amine();
             // check this atom is a guanidino C (sp2 valence, C.cat)
             bool is_guanidino();
             //check this carbon is a amidine C (sp2 valence)
             bool is_amidine();

             // check if the hetero atom belongs to aromatic ring 
             bool is_aromatic_hetero();
             // check this atom is planar nitrogen
             bool is_planar_amide();
                                           
             // Now define some type of H-bond acceptor atoms with different hybridization
             // 1. Acceptor
             // check this atom is a sp2 oxygen(or sulfur) of H-bond acceptor
             //(-CO-, -SO-, -SO2-, -PO2-, -PO3, -CS-)
             inline bool is_sp2_oxygen_acceptor(){return (symbol_type_ == "O.2" || symbol_type_ == "S.2") && num_of_bonded_atom_ == 1;}
             // check this atom is a sp3 oxygen(or sulfur) of HB acceptor (caution: this type of acceptor is weak so we 
			 // discard those in localized systems (e.g. oxygen bridge 2 pi atoms)
             bool is_sp3_oxygen_acceptor();
             // check this atom is sp2 planar nitrogen of H-bond acceptor
             // (=N-, aromatic nitrogen)
             inline bool is_sp2_nitrogen_acceptor(){return ((symbol_type_ == "N.2" || symbol_type_ == "N.ar") && num_of_bonded_atom_ == 2);}
             // check this atom is sp nitrogen of H-bond acceptor
             // (-C-=N)
             inline bool is_sp_nitrogen_acceptor(){return symbol_type_ == "N.1" && num_of_bonded_atom_ == 1;}
             // check this atom is a carboxyl oxygen of H-bond acceptor
             // -COO
             inline bool is_carboxyl_acceptor(){return is_carboxyl_oxygen();}
             
             // 2. Donor
             // check if this atom is sp3 nitrogen of H-bond donor
             // -NH2, -NH3+, -NH-
             inline bool is_sp3_nitrogen_donor(){return bond_to_hydrogen() && (symbol_type_ == "N.3" || symbol_type_ == "N.4");}
             // check if this atom is sp2 nitrogen of H-bond donor
             // -CONH-, triangle planar nitrogen, =NH,
             inline bool is_sp2_nitrogen_donor(){return (bond_to_hydrogen()&& (symbol_type_ == "N.2" || symbol_type_ == "N.am"
                                                 || symbol_type_ == "N.ar" || symbol_type_ == "N.pl3"));}
             // check if this atom is a sp3 oxygen(sulfur) of H-bond donor
             // -OH, -SH
             inline bool is_sp3_oxygen_donor(){return is_hydroxyl() || is_sulfydryl();}
             
             // check if this atom is positively charged (not ionizable)
             // generally speaking, nitrogen or phosphorus bearing valence 4
             bool is_pos_charged();
             // check if this atom is negatively charged (not ionizable)
             // generally speaking, sp3 hybrition oxygen or sulfur bearing bond order 1
             bool is_neg_charged();
             
             // rotational mark related methods
             void set_rotation_mark(ATOM* anchor);
             inline bool get_rotation_mark()
			 {
				 return rotation_mark_;
			 }
             inline void clear_rotation_mark(){rotation_mark_ = false;}
             
             //reset the postion of the atom to its original value
             //inline void reset_position(){_p = orig_position_;}
			 //judging whether the two atoms are in the same ring or not  


            //! \return The sum of the bond orders of the bonds to the atom (i.e. double bond = 2...)
			 unsigned int  BOSum();
			 int get_atomic_num() const;
			
			
      private:
              int id;
			  vector3 _p;
              vector3 orig_position_;
              vector3 _f;
              string name;
              string symbol_type_;
              string element_;
              int type;
			  int mmff94_type;
			  string mmff94_symbol_type_;
              int subst_id;
              string subst_name;
              float charge;
			  float mmff94_partial_charge;
			  float radius_;
			  float weight_;
              //int valence_; // 1-sp, 2-sp2, 3-sp3, 4-delocalized-pi
              int num_of_bond_;
			  int num_of_rotor_bond_;
              int num_of_bonded_atom_;
              int num_of_bonded_heavy_atom_;
              // vector holding neighbor bonds
              BONDVec bond_;
              // vector holding neighbor bonded atoms
              ATOMVec atom_;
              // a mark flag to indicate this atom is part of the rotational fragment
              bool rotation_mark_;
			  // a mark flag to indicate this atom is an aromatic ring atom
			  bool AromaticRingAtom_;
};

#endif        




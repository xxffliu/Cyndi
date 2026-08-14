#ifndef FORCEFIELD_H
#define FORCEFIELD_H

#ifndef MOL_H
#include "Mol.h"
#endif

#ifndef FFPARAMETER_H
#include "FFParameter.h"
#endif

#include "utility.h"

#include <vector>

// avogadro constant
static const double AVOGADRO = 6.0221367e+23;
// Joule per calorie
static const double JOULE_PER_CAL = 4.1868;
// conversion from Kcal/(mol A) into Newton
static const double FORCE_FACTOR = JOULE_PER_CAL * 1000 * 1E10 / AVOGADRO;

class FFComponent;
//class FRAGMENT;
static const int DEFAULT_UPDATEFREQUENCY = 10;
// Avogadro constant.
//const double	AVOGADRO        	= 6.0221367E+23L;   	 // 1 / mol
class ForceField{
      public:
             friend class FFComponent;
             
             // constructors and destructor
             ForceField();
             ForceField(MOL& mol);
             ForceField(const ForceField& ff);
             virtual ~ForceField();
             virtual void clear();
             // assignment operation
             ForceField& operator=(const ForceField& ff);
             // debug
             bool isValid() const;
             
             // setup methods: set up the forceld and its components
             bool setup(MOL& mol_);
             // Force Field specific setup, this method is called by setup
             virtual bool specific_setup();
             // return the list of unassigned atoms
             ATOMVec& get_unassigned_atoms();
             // add atoms with unassigned force field types
             void add_unassigned_atom(ATOM*);
             // set and get the force field name
             void set_ff_name(const string& name);
             string get_ff_name();
             
             // return the number of atoms stored in force field
             int get_num_atoms();
             // return a (const) pointer to the atom vector
             inline ATOMVec& get_atoms(){return atom_vec_;}
             inline const ATOMVec& get_atoms() const{ return atom_vec_;}
             // return a (const) pointer to the bond vector
             inline BONDVec& get_bonds(){return bond_vec_;}
             inline const BONDVec& get_bonds() const{ return bond_vec_;}
             // return a (const) pointer to the molecule
             inline MOL* get_mol(){ return mol_;}
             inline const MOL* get_mol() const {return mol_;}
             
             //return a pointer to the force field parameter file;
             
             FFParameter& get_parameters();
             //insert a new force field component into the force field component list
             void insert_component(FFComponent* ffp);
             // remove a component from the force field component list
             void remove_component(const FFComponent* ffp);
             void remove_component(const string& component);
             // return a pointer to a specific force field component
             FFComponent* get_component(const int id) const;
             FFComponent* get_component(const string& component) const;
             
             // return the sum of enenrgies registered in FFComponent in kj/mol
             double get_energy() const;
             
             // update the current enegy
             double update_energy();
             
             // update the sum of forces imposed on each atoms by every ff component
             void update_forces();
             
             // calculate rms of current gradient
             double get_rms_gradient() const;
             
             // return the update fequency for pair lists
             virtual int get_update_frequency() const;
             
             // update internal data structure
             virtual void update();
			 //judging whether the atom a and the atom b are in the same ring
			 bool  In_the_sameRing(ATOM* atom1,ATOM* atom2);
			 
      protected:
                MOL* mol_;
				FFParameter parameter_;
                ATOMVec atom_vec_;
                BONDVec bond_vec_;
                bool valid_;
                string ff_name_;
                double energy_;
                vector<FFComponent*> component_;
                ATOMVec unassignedAtoms_;
				
};
#endif
                
             
             
             
             
             
                
                
